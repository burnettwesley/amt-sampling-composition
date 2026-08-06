# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 56_synthetic_uniform_hours.R
# ============================
# Outcome-free rerun of the scripts/45 synthetic-outcome experiment (baseline
# 2 km / +/-72 h cell) and the scripts/55 parameter sweep, for both cities,
# changing exactly ONE thing: each synthetic site's hour is drawn UNIFORMLY
# from the set of distinct UTC clock hours in which the city's fleet recorded
# at least one valid NOx reading ("observed collection hours"), instead of
# from the Oakland crash hour-of-day distribution used by scripts/45/55.
# This script reads NO outcome data of any kind (no crashes, no violations),
# so the AMT paper's bright line holds within the experiment itself.
#
# Site locations, collection days, and within-hour second offsets are
# bit-identical to the archived scripts/45 run: seed 20260720, SF processed
# first, and scripts/45's RNG call sequence is replicated exactly, including
# a discarded burn draw (sample.int(1217L, N_SYN, replace = TRUE)) that
# consumes the same generator output as scripts/45's draw from its
# 1,217-element crash-hour vector (sample(x, ...) delegates to
# sample.int(length(x), ...), so consumption depends only on the length).
# Only the integer 1217 -- a row count -- is retained from the old design; no
# outcome file is read. The NEW uniform hour draws come from isolated
# sub-streams (seed 20260721 for SF, 20260722 for Oakland) with the main
# stream saved and restored around them, so the main stream -- and hence both
# cities' site draws -- match scripts/45 exactly. Only the times change.
#
# Matching (blind and class-constrained) uses scripts/55's fast superset
# implementation, which reproduced the archived scripts/45 baseline exactly
# (see the REPRODUCED checks in output/tables/amt_sensitivity.txt).
#
# OUTPUT: output/tables/amt_synthetic_uniform.txt
#         (sf_siting.txt and amt_sensitivity.txt are NOT touched)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(lubridate)
})
if (!file.exists("scripts/45_sf_siting.R"))
  stop("Run this script from the repository root.")

set.seed(20260720)   # same seed as scripts/45; SF is processed first, as
                     # there, so all site/day/second draws match that run
HOUR_SEED <- c(SF = 20260721, Oakland = 20260722)

UTM10N  <- 32610
N_SYN   <- 2000
N_BURN  <- 1217L     # length of scripts/45's crash-hour vector (burn only)
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

# -- loading (verbatim from scripts/45/55) ------------------------------------
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
# (st_sample; sample(days); BURN in place of the old crash-hour draw; runif)
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
  burn <- sample.int(N_BURN, N_SYN, replace = TRUE)   # discarded; see header
  s_sec <- runif(N_SYN, 0, 3600)
  list(sxy = sxy, s_class = s_class, s_date = s_date, s_sec = s_sec,
       mp = pts %>% select(ts, x, y, road_class, nox))
}

# -- parameter sweep: verbatim from scripts/55, plus paired means -------------
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
           mn_b    = mean(val_b[both, cc]),
           mn_c    = mean(val_c[both, cc]),
           gap     = mean(val_b[both, cc]) - mean(val_c[both, cc]),
           n_both  = sum(both))
  })
}

# -- run both cities (SF FIRST -- do not reorder: RNG stream must match 45) ---
res <- list(); hrs <- list()
for (city in c("SF", "Oakland")) {
  t_start <- Sys.time()
  cat(sprintf("[%s] loading readings ...\n", city))
  d <- load_city(city)
  cat(sprintf("[%s] %s readings; generating synthetic sites ...\n",
              city, format(nrow(d$pts), big.mark = ",")))
  sites <- gen_sites(d)
  # NEW: hours uniform over the city's observed collection hours, drawn on an
  # isolated sub-stream so the main stream stays identical to scripts/45.
  obs_hours <- sort(unique(hour(sites$mp$ts)))
  saved <- .Random.seed
  set.seed(HOUR_SEED[[city]])
  s_hour <- sample(obs_hours, N_SYN, replace = TRUE)
  .Random.seed <- saved
  sites$s_ts <- as.POSIXct(paste(sites$s_date), tz = "UTC") +
    s_hour * 3600 + sites$s_sec
  hrs[[city]] <- obs_hours
  cat(sprintf("[%s] %d observed collection hours (UTC): %s\n",
              city, length(obs_hours), paste(obs_hours, collapse = " ")))
  cat(sprintf("[%s] parameter sweep (%d cells x blind/constrained) ...\n",
              city, nrow(CELLS)))
  res[[city]] <- run_sweep(sites)
  cat(sprintf("[%s] done in %.1f min\n", city,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  rm(d, sites); invisible(gc())
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

out <- file("output/tables/amt_synthetic_uniform.txt", open = "wt")
w <- function(...) cat(..., "\n", sep = "", file = out)
w("AMT synthetic experiment, uniform collection-hour draws -- scripts/56 -- ",
  format(Sys.Date()))
w("Synthetic sites: 2000 per city, non-freeway segments, length-weighted;")
w("site locations, collection days, and within-hour offsets identical to the")
w("archived scripts/45 run (seed 20260720, SF first; RNG stream alignment")
w("via a discarded burn draw -- see the script header). Hours drawn")
w("UNIFORMLY from the set of distinct UTC clock hours with at least one")
w("valid NOx reading in that city. NO outcome data of any kind are read:")
w("no crashes, no violations. Matching rule identical to scripts/45/55")
w("(combined metric d = sqrt((m/500)^2 + (hr/6)^2); hard caps as labeled;")
w("'*' marks the baseline cell, 2 km / +/-72 h).")
w("Compo-infl: blind-assigned class mix minus own-class mix, valued at")
w("citywide class means of log(NOx) (log points). Paired gap: blind minus")
w("constrained mean log(NOx) on the both-matched subset.")
w("")
for (city in c("SF", "Oakland")) {
  b <- res[[city]] %>% filter(R == 2000, W == 72)
  w(sprintf("== %s: synthetic-outcome matching experiment (N = %d non-freeway sites) ==",
            city, N_SYN))
  w(sprintf("  observed collection hours (UTC, n = %d): %s",
            length(hrs[[city]]), paste(hrs[[city]], collapse = " ")))
  w(sprintf("  match rate: blind %.1f%%, class-constrained %.1f%%",
            b$match_b, b$match_c))
  w(sprintf("  blind: share assigned a FREEWAY reading: %.1f%%", b$fw))
  w(sprintf("  blind: share assigned a reading of the site's own class: %.1f%%",
            b$own))
  w(sprintf("  composition-implied inflation (blind class mix vs own-class mix, at citywide class means): %.2f log points",
            b$infl))
  w(sprintf("  paired blind-vs-constrained gap (both-matched subset, n = %d): blind %.2f, constrained %.2f (gap %.2f log points)",
            b$n_both, b$mn_b, b$mn_c, b$gap))
  w("")
  w("  parameter sweep:")
  for (ln in fmt_sweep(res[[city]])) w(ln)
  w("")
}
close(out)
cat(readLines("output/tables/amt_synthetic_uniform.txt"), sep = "\n")
