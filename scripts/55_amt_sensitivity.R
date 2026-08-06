# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 55_amt_sensitivity.R
# ====================
# Sensitivity sweeps for the AMT companion paper, extending scripts/44-45.
# Three panels per city (Oakland, San Francisco), readings + synthetic sites
# only -- NO outcome data (no crashes, no violations) per the AMT bright line.
#
# Panel 1 -- SYNTHETIC-EXPERIMENT PARAMETER SWEEP: rerun the scripts/45
#   blind-vs-class-constrained matching experiment varying the spatial radius
#   (1 km, 2 km baseline, 3 km) and the temporal window (+/-24 h, +/-72 h
#   baseline, +/-168 h), holding the other cap at baseline. The combined
#   metric d = sqrt((m/500)^2 + (hr/6)^2) is held fixed; only the hard caps
#   vary. Seed and city order (SF first) are identical to scripts/45, so the
#   RNG stream -- and hence the 2,000 synthetic non-freeway sites per city --
#   are identical to the archived baseline, and identical across sweep cells.
#   Tie-breaking replicates scripts/45's which.min in original row order.
#
# Panel 2 -- CELL-SIZE ROBUSTNESS: median revisit days per cell by road class
#   at 50 m, 100 m (baseline), and 200 m grid cells.
#
# Panel 3 -- REVISIT DISTRIBUTIONS: p25/p50/p75/p90 of distinct visit days
#   per 100 m cell by road class; compact two-panel figure.
#
# OUTPUT: output/tables/amt_sensitivity.txt
#         output/figures/fig_amt_revisit_dist.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(lubridate)
})
if (!file.exists("scripts/45_sf_siting.R"))
  stop("Run this script from the repository root.")

set.seed(20260720)   # same seed as scripts/45; SF is processed first, as
                     # there, so all RNG draws match the archived run

UTM10N  <- 32610
N_SYN   <- 2000
CA_FILE <- Sys.getenv("AQ_CA_FILE", "data/raw/aclima/California_201605_201709_GoogleAclimaAQ.txt")
SFBB    <- c(lat_min = 37.708, lat_max = 37.833, lon_min = -122.515, lon_max = -122.355)
cls_lab <- c("local", "collector", "arterial", "freeway")

CELLS <- tibble(
  label = c("1 km / +/-72 h",
            "2 km / +/-72 h *",
            "3 km / +/-72 h",
            "2 km / +/-24 h",
            "2 km / +/-168 h"),
  R = c(1000, 2000, 3000, 2000, 2000),
  W = c(  72,   72,   72,   24,  168))
R_MAX <- max(CELLS$R); W_MAX <- max(CELLS$W)

SIZES <- tibble(size = c("50 m", "100 m", "200 m"),
                lon  = c(0.0006, 0.0012, 0.0024),
                lat  = c(0.00045, 0.0009, 0.0018))

# archived scripts/45 baseline (output/tables/sf_siting.txt, 2026-07-20):
# blind match %, freeway-assigned %, own-class %, composition inflation
ARCHIVE <- list(Oakland = c(44.5, 64.5, 11.7, 1.05),
                SF      = c(42.0, 36.7, 21.4, 0.81))
ARCH_MED <- list(  # 100 m medians: arterial, collector, freeway, local
  Oakland = c(arterial = 33, collector = 26, freeway = 65, local = 16),
  SF      = c(arterial =  9, collector =  5, freeway = 34, local =  3))

# Oakland crash hour-of-day distribution (UTC), for outcome-time draws only
# (design input inherited from scripts/45; no outcome variables are used)
cr_hours <- read_csv("data/processed/nearest_matched.csv",
                     show_col_types = FALSE) %>%
  mutate(h = hour(as.POSIXct(timestamp, tz = "UTC"))) %>% pull(h)

# -- loading (verbatim from scripts/45) ---------------------------------------
load_city <- function(which) {
  if (which == "SF") {
    pts <- read_csv(CA_FILE,
                    col_select = c(Date_Time, Latitude, Longitude, NO2, NO),
                    show_col_types = FALSE) %>%
      mutate(ts = ymd_hms(Date_Time, tz = "UTC"),
             NO2 = suppressWarnings(as.numeric(NO2)),
             NO  = suppressWarnings(as.numeric(NO))) %>%
      filter(!is.na(Latitude), !is.na(Longitude), !is.na(ts),
             Latitude >= SFBB["lat_min"], Latitude <= SFBB["lat_max"],
             Longitude >= SFBB["lon_min"], Longitude <= SFBB["lon_max"])
    rc <- read_csv("data/processed/sf_reading_road_class.csv.gz",
                   show_col_types = FALSE)
    stopifnot(nrow(rc) == nrow(pts))
    pts$road_class <- rc$road_class
    roads <- st_read("data/processed/sf_osm_roads.gpkg", quiet = TRUE)
  } else {
    pts <- read_csv("data/raw/aclima/Oakland_201505_201605_GoogleAclimaAQ.txt",
                    col_select = c(Date_Time, Latitude, Longitude, NO2, NO),
                    show_col_types = FALSE) %>%
      mutate(ts = ymd_hms(Date_Time, tz = "UTC"),
             NO2 = suppressWarnings(as.numeric(NO2)),
             NO  = suppressWarnings(as.numeric(NO))) %>%
      filter(!is.na(Latitude), !is.na(Longitude), !is.na(ts))
    rc <- read_csv("data/processed/reading_road_class.csv.gz",
                   show_col_types = FALSE)
    stopifnot(nrow(rc) == nrow(pts))
    pts$road_class <- rc$road_class
    roads <- st_read("data/processed/osm_roads.gpkg", quiet = TRUE)
  }
  xy <- pts %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(UTM10N) %>% st_coordinates()
  pts$x <- xy[, 1]; pts$y <- xy[, 2]
  pts <- pts %>% mutate(nox = NO2 + NO, date = as.Date(ts))
  list(pts = pts, roads = roads)
}

# -- synthetic sites: RNG call sequence verbatim from scripts/45 --------------
gen_sites <- function(d) {
  pts <- d$pts %>% filter(!is.na(nox), nox > 0)
  roads <- d$roads
  nonf <- roads %>% filter(road_class != "freeway")
  samp <- st_sample(nonf, size = ceiling(N_SYN * 1.3), type = "random")
  samp <- st_cast(samp, "POINT")
  sxy_full <- st_coordinates(samp)[, 1:2, drop = FALSE]
  good <- which(is.finite(sxy_full[, 1]) & is.finite(sxy_full[, 2]))
  stopifnot(length(good) >= N_SYN)
  good <- good[seq_len(N_SYN)]
  samp <- samp[good]
  sxy <- sxy_full[good, , drop = FALSE]
  sidx <- st_nearest_feature(st_sf(geometry = samp), roads)
  s_class <- roads$road_class[sidx]
  stopifnot(length(s_class) == N_SYN, nrow(sxy) == N_SYN)
  days <- sort(unique(pts$date))
  s_date <- sample(days, N_SYN, replace = TRUE)
  s_hour <- sample(cr_hours, N_SYN, replace = TRUE)
  s_ts <- as.POSIXct(paste(s_date), tz = "UTC") + s_hour * 3600 +
    runif(N_SYN, 0, 3600)
  list(sxy = sxy, s_class = s_class, s_ts = s_ts,
       mp = pts %>% select(ts, x, y, road_class, nox))
}

# -- Panel 1: parameter sweep -------------------------------------------------
# One superset extraction per site (<= R_MAX, +/-W_MAX); each sweep cell is a
# subset. Candidate arithmetic (sp, tmp, d) is identical to scripts/45;
# among tied minima the candidate earliest in original row order is chosen,
# replicating scripts/45's which.min over an original-order candidate frame.
run_sweep <- function(sites) {
  mp <- sites$mp
  orig <- seq_len(nrow(mp))
  o <- order(as.numeric(mp$ts))
  tsn <- as.numeric(mp$ts)[o]
  px <- mp$x[o]; py <- mp$y[o]
  cint <- match(mp$road_class, cls_lab)[o]
  lnx <- log(mp$nox)[o]
  orig <- orig[o]
  s_cint <- match(sites$s_class, cls_lab)
  t0v <- as.numeric(sites$s_ts)
  nC <- nrow(CELLS)
  cls_b <- matrix(NA_integer_, N_SYN, nC); val_b <- matrix(NA_real_, N_SYN, nC)
  cls_c <- cls_b; val_c <- val_b
  for (i in seq_len(N_SYN)) {
    t0 <- t0v[i]
    lo <- findInterval(t0 - W_MAX * 3600, tsn, left.open = TRUE) + 1L
    hi <- findInterval(t0 + W_MAX * 3600, tsn)
    if (hi < lo) next
    ii <- lo:hi
    dx <- px[ii] - sites$sxy[i, 1]; dy <- py[ii] - sites$sxy[i, 2]
    sp <- sqrt(dx * dx + dy * dy)
    k <- sp <= R_MAX
    if (!any(k)) next
    spk <- sp[k]
    tk  <- abs(tsn[ii][k] - t0) / 3600
    ck  <- cint[ii][k]; lk <- lnx[ii][k]; ok <- orig[ii][k]
    for (cc in seq_len(nC)) {
      sel <- which(spk <= CELLS$R[cc] & tk <= CELLS$W[cc])
      if (!length(sel)) next
      dd <- sqrt((spk[sel] / 500)^2 + (tk[sel] / 6)^2)
      tie <- sel[dd == min(dd)]
      j <- tie[which.min(ok[tie])]
      cls_b[i, cc] <- ck[j]; val_b[i, cc] <- lk[j]
      selc <- sel[ck[sel] == s_cint[i]]
      if (length(selc)) {
        ddc <- sqrt((spk[selc] / 500)^2 + (tk[selc] / 6)^2)
        tiec <- selc[ddc == min(ddc)]
        jc <- tiec[which.min(ok[tiec])]
        cls_c[i, cc] <- ck[jc]; val_c[i, cc] <- lk[jc]
      }
    }
  }
  cm <- mp %>% group_by(road_class) %>%
    summarise(m = mean(log(nox)), .groups = "drop")
  mvec <- setNames(cm$m, cm$road_class)
  map_dfr(seq_len(nC), function(cc) {
    okb <- !is.na(cls_b[, cc]); okc <- !is.na(cls_c[, cc])
    both <- okb & okc
    tibble(label = CELLS$label[cc], R = CELLS$R[cc], W = CELLS$W[cc],
           match_b = 100 * mean(okb),
           fw      = 100 * mean(cls_b[okb, cc] == 4L),
           own     = 100 * mean(cls_b[okb, cc] == s_cint[okb]),
           infl    = mean(mvec[cls_lab[cls_b[okb, cc]]]) -
                     mean(mvec[sites$s_class[okb]]),
           match_c = 100 * mean(okc),
           gap     = mean(val_b[both, cc]) - mean(val_c[both, cc]),
           n_both  = sum(both))
  })
}

# -- Panels 2-3: revisit frequency by cell size -------------------------------
# Cell definition and modal-class rule as in scripts/45 (100 m baseline:
# floor(lon/0.0012), floor(lat/0.0009); ties resolved alphabetically, matching
# scripts/45's stable sort of the class table). Uses ALL readings (incl. NA
# NOx), as in scripts/45's diagnostics.
revisit_by_size <- function(pts, lon_div, lat_div) {
  p <- pts %>% transmute(
    cell = paste(floor(Longitude / lon_div), floor(Latitude / lat_div)),
    road_class, date)
  modal <- p %>% count(cell, road_class) %>%
    arrange(cell, desc(n), road_class) %>%
    distinct(cell, .keep_all = TRUE) %>%
    select(cell, class = road_class)
  nd <- p %>% distinct(cell, date) %>% count(cell, name = "n_days")
  left_join(nd, modal, by = "cell")
}

# -- run both cities (SF FIRST -- do not reorder: RNG stream must match 45) ---
res <- list(); rv_fig <- list()
for (city in c("SF", "Oakland")) {
  t_start <- Sys.time()
  cat(sprintf("[%s] loading readings ...\n", city))
  d <- load_city(city)
  cat(sprintf("[%s] %s readings; generating synthetic sites ...\n",
              city, format(nrow(d$pts), big.mark = ",")))
  sites <- gen_sites(d)
  cat(sprintf("[%s] parameter sweep (%d cells x blind/constrained) ...\n",
              city, nrow(CELLS)))
  sw <- run_sweep(sites)
  cat(sprintf("[%s] revisit panels ...\n", city))
  rv <- map(seq_len(nrow(SIZES)),
            \(s) revisit_by_size(d$pts, SIZES$lon[s], SIZES$lat[s]))
  names(rv) <- SIZES$size
  med_tab <- imap_dfr(rv, ~ .x %>% group_by(class) %>%
                        summarise(med = median(n_days),
                                  n_cells = n(), .groups = "drop") %>%
                        mutate(size = .y))
  q_tab <- rv[["100 m"]] %>% group_by(class) %>%
    summarise(n_cells = n(),
              p25 = quantile(n_days, .25), p50 = quantile(n_days, .50),
              p75 = quantile(n_days, .75), p90 = quantile(n_days, .90),
              .groups = "drop")
  rv_fig[[city]] <- rv[["100 m"]] %>% select(class, n_days) %>%
    mutate(city = city)
  res[[city]] <- list(sweep = sw, med = med_tab, q = q_tab)
  cat(sprintf("[%s] done in %.1f min\n", city,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  rm(d, sites, rv); invisible(gc())
}

# -- write table --------------------------------------------------------------
fmt_sweep <- function(sw) {
  hdr <- sprintf("  %-18s %11s %13s %9s %10s | %12s %14s",
                 "cell", "blind-match", "freeway-share", "own-class",
                 "compo-infl", "constr-match", "paired-gap (n)")
  rows <- sprintf("  %-18s %10.1f%% %12.1f%% %8.1f%% %10.2f | %11.1f%% %9.2f (%d)",
                  sw$label, sw$match_b, sw$fw, sw$own, sw$infl,
                  sw$match_c, sw$gap, sw$n_both)
  c(hdr, rows)
}
fmt_med <- function(med) {
  w <- med %>% select(size, class, med) %>%
    pivot_wider(names_from = class, values_from = med)
  hdr <- sprintf("  %-7s %8s %10s %9s %8s",
                 "size", "local", "collector", "arterial", "freeway")
  rows <- sprintf("  %-7s %8g %10g %9g %8g",
                  w$size, w$local, w$collector, w$arterial, w$freeway)
  c(hdr, rows)
}
fmt_q <- function(q) {
  q <- q %>% mutate(class = factor(class, levels = cls_lab)) %>% arrange(class)
  hdr <- sprintf("  %-10s %8s %6s %6s %6s %6s",
                 "class", "n_cells", "p25", "p50", "p75", "p90")
  rows <- sprintf("  %-10s %8d %6.1f %6.1f %6.1f %6.1f",
                  as.character(q$class), q$n_cells, q$p25, q$p50, q$p75, q$p90)
  c(hdr, rows)
}

base_check <- function(city) {
  b <- res[[city]]$sweep %>% filter(R == 2000, W == 72)
  got <- c(round(b$match_b, 1), round(b$fw, 1), round(b$own, 1),
           round(b$infl, 2))
  tgt <- ARCHIVE[[city]]
  ok <- isTRUE(all(abs(got - tgt) < 1e-9))
  list(ok = ok,
       msg = sprintf(
         "  baseline (2 km / +/-72 h) vs archived sf_siting.txt: %s (got %.1f%% / %.1f%% / %.1f%% / %.2f; archived %.1f%% / %.1f%% / %.1f%% / %.2f)",
         if (ok) "REPRODUCED" else "MISMATCH",
         got[1], got[2], got[3], got[4], tgt[1], tgt[2], tgt[3], tgt[4]))
}
med_check <- function(city) {
  m <- res[[city]]$med %>% filter(size == "100 m")
  got <- setNames(m$med, m$class)[names(ARCH_MED[[city]])]
  ok <- isTRUE(all(got == ARCH_MED[[city]]))
  list(ok = ok,
       msg = sprintf("  100 m medians vs archived sf_siting.txt: %s",
                     if (ok) "REPRODUCED" else "MISMATCH"))
}

out <- file("output/tables/amt_sensitivity.txt", open = "wt")
w <- function(...) cat(..., "\n", sep = "", file = out)
w("AMT sensitivity sweeps -- scripts/55 -- ", format(Sys.Date()))
w("Seed 20260720, SF processed first (both as in scripts/45), so the RNG")
w("stream and hence the 2,000 synthetic non-freeway sites per city are")
w("identical to the archived baseline and identical across sweep cells.")
w("Panel 1 varies only the hard caps (spatial radius, temporal window);")
w("the combined metric d = sqrt((m/500)^2 + (hr/6)^2) is held fixed.")
w("'*' marks the scripts/45 baseline cell (2 km / +/-72 h).")
w("Compo-infl: blind-assigned class mix minus own-class mix, valued at")
w("citywide class means of log(NOx) (log points). Paired gap: blind minus")
w("constrained mean log(NOx) on the both-matched subset. No outcome data.")
w("")
for (city in c("Oakland", "SF")) {
  bc <- base_check(city); mc <- med_check(city)
  w(sprintf("==================== %s ====================", city))
  w("")
  w("Panel 1: synthetic-experiment parameter sweep (N = 2000 non-freeway sites)")
  for (ln in fmt_sweep(res[[city]]$sweep)) w(ln)
  w(bc$msg)
  w("")
  w("Panel 2: median revisit days per cell by road class, by cell size")
  for (ln in fmt_med(res[[city]]$med)) w(ln)
  w(mc$msg)
  w("")
  w("Panel 3: distinct visit days per 100 m cell -- quartiles by road class")
  for (ln in fmt_q(res[[city]]$q)) w(ln)
  w("")
}
close(out)
cat(readLines("output/tables/amt_sensitivity.txt"), sep = "\n")

# -- figure -------------------------------------------------------------------
figd <- bind_rows(rv_fig) %>%
  mutate(class = factor(class, levels = cls_lab),
         city = factor(city, levels = c("Oakland", "SF"),
                       labels = c("Oakland (Jun 2015 - May 2016)",
                                  "San Francisco (Jun 2016 - Sep 2017)")))
p <- ggplot(figd, aes(class, n_days)) +
  geom_boxplot(outlier.size = 0.25, outlier.alpha = 0.25,
               fill = "grey88", linewidth = 0.3) +
  scale_y_log10() +
  facet_wrap(~ city, nrow = 1) +
  labs(x = NULL, y = "Distinct visit days per 100 m cell") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave("output/figures/fig_amt_revisit_dist.png", p,
       width = 6.5, height = 2.6, dpi = 300, bg = "white")
cat("Saved output/figures/fig_amt_revisit_dist.png\n")
