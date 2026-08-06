# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 65_ufp_replication.R
# ====================
# Registered-prediction test: ultrafine particles (particle number, PN) are a
# PRIMARY traffic pollutant, emitted directly by combustion, so the
# siting-composition mechanism predicts a POSITIVE artifact with the same
# sign as NOx -- PN highest on freeways, blind joins OVERSTATING PN at
# non-freeway sites -- in contrast to the sign-reversed O3 result
# (scripts/64).
#
# DESIGN: SF only (the California_201605_201709 release carries PN1-PN5
# columns; the Oakland release does not). We replicate the scripts/56 SF
# synthetic experiment BIT-IDENTICALLY -- same seed (20260720, SF processed
# first, so the RNG stream needs no Oakland pass), same 2,000 non-freeway
# sites, same collection days, same uniform hour draws from the isolated
# sub-stream (seed 20260721), same discarded burn draw, same matching pool
# (readings with valid NOx > 0) and the same combined-metric join. The JOIN
# IS UNCHANGED: each synthetic site is matched to exactly the reading it
# received in scripts/56; we then read off that reading's PN instead of its
# NOx. A hard reproduction check confirms the baseline-cell NOx numbers
# equal the archived scripts/56 values before any PN result is written.
#
# PN MEASURE AND SUPPORT (documented in output): the release reports
# particle number in five size bins (PN1-PN5, particle counts; PN1 is the
# smallest/dominant bin). Coverage in the SF box is 91.4% for all five bins
# jointly (91.4% for PN1 alone; the four coarser bins add no coverage cost),
# so we use TOTAL PN = PN1 + ... + PN5 as the measured pollutant; requiring
# all five bins loses essentially nothing relative to PN1 alone, and the
# total is the natural particle-number concentration. Non-positive totals
# are essentially absent (0.00% of valid readings), so the support permits
# logs; log(PN) is reported as the primary scale for direct comparability
# with the NOx log-point figures, with levels alongside.
#
# NO outcome data of any kind are read (no crashes, no violations).
#
# OUTPUT: output/tables/ufp_replication.txt
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(lubridate)
})
if (!file.exists("scripts/45_sf_siting.R"))
  stop("Run this script from the repository root.")

set.seed(20260720)   # same seed as scripts/45/56; SF is processed first
HOUR_SEED_SF <- 20260721

UTM10N  <- 32610
N_SYN   <- 2000
N_BURN  <- 1217L     # length of scripts/45's crash-hour vector (burn only)
CA_FILE <- Sys.getenv("AQ_CA_FILE", "data/raw/aclima/California_201605_201709_GoogleAclimaAQ.txt")
SFBB    <- c(lat_min = 37.708, lat_max = 37.833, lon_min = -122.515, lon_max = -122.355)
cls_lab <- c("local", "collector", "arterial", "freeway")

# archived scripts/56 SF baseline (2 km / +/-72 h), for the reproduction check
ARCH <- list(match_b = 41.9, match_c = 34.3, fw = 36.0, own = 21.6,
             infl = 0.81, mn_b = 2.71, mn_c = 2.56, gap = 0.15, n_both = 686L)

CELLS <- tibble(
  label = c("1 km / +/-72 h",
            "2 km / +/-72 h *",
            "3 km / +/-72 h",
            "2 km / +/-24 h",
            "2 km / +/-168 h"),
  R = c(1000, 2000, 3000, 2000, 2000),
  W = c(  72,   72,   72,   24,  168))
R_MAX <- max(CELLS$R); W_MAX <- max(CELLS$W)

# -- loading (scripts/56 verbatim, plus the PN columns) -----------------------
load_sf <- function() {
  pts <- read_csv(CA_FILE,
                  col_select = c(Date_Time, Latitude, Longitude,
                                 PN1, PN2, PN3, PN4, PN5, NO2, NO),
                  show_col_types = FALSE) %>%
    mutate(ts = ymd_hms(Date_Time, tz = "UTC"),
           across(c(PN1, PN2, PN3, PN4, PN5, NO2, NO),
                  ~ suppressWarnings(as.numeric(.x)))) %>%
    filter(!is.na(Latitude), !is.na(Longitude), !is.na(ts),
           Latitude >= SFBB["lat_min"], Latitude <= SFBB["lat_max"],
           Longitude >= SFBB["lon_min"], Longitude <= SFBB["lon_max"])
  rc <- read_csv("data/processed/sf_reading_road_class.csv.gz",
                 show_col_types = FALSE)
  stopifnot(nrow(rc) == nrow(pts))
  pts$road_class <- rc$road_class
  roads <- st_read("data/processed/sf_osm_roads.gpkg", quiet = TRUE)
  xy <- pts %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(UTM10N) %>% st_coordinates()
  pts$x <- xy[, 1]; pts$y <- xy[, 2]
  pts <- pts %>% mutate(nox = NO2 + NO,
                        PN  = PN1 + PN2 + PN3 + PN4 + PN5,  # NA if any bin NA
                        date = as.Date(ts))
  list(pts = pts, roads = roads)
}

# -- synthetic sites: RNG call sequence verbatim from scripts/45/56 -----------
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
       mp = pts %>% select(ts, x, y, road_class, nox, PN))
}

# -- sweep: scripts/56 join verbatim, but store matched-row INDICES -----------
run_sweep_idx <- function(sites) {
  mp <- sites$mp
  orig <- seq_len(nrow(mp))
  o <- order(as.numeric(mp$ts))
  tsn <- as.numeric(mp$ts)[o]
  px <- mp$x[o]; py <- mp$y[o]
  cint <- match(mp$road_class, cls_lab)[o]
  orig <- orig[o]
  s_cint <- match(sites$s_class, cls_lab)
  t0v <- as.numeric(sites$s_ts)
  nC <- nrow(CELLS)
  idx_b <- matrix(NA_integer_, N_SYN, nC); idx_c <- idx_b
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
    ck  <- cint[ii][k]; ok <- orig[ii][k]
    for (cc in seq_len(nC)) {
      sel <- which(spk <= CELLS$R[cc] & tk <= CELLS$W[cc])
      if (!length(sel)) next
      dd <- sqrt((spk[sel] / 500)^2 + (tk[sel] / 6)^2)
      tie <- sel[dd == min(dd)]
      j <- tie[which.min(ok[tie])]
      idx_b[i, cc] <- ok[j]
      selc <- sel[ck[sel] == s_cint[i]]
      if (length(selc)) {
        ddc <- sqrt((spk[selc] / 500)^2 + (tk[selc] / 6)^2)
        tiec <- selc[ddc == min(ddc)]
        jc <- tiec[which.min(ok[tiec])]
        idx_c[i, cc] <- ok[jc]
      }
    }
  }
  list(idx_b = idx_b, idx_c = idx_c)
}

# =============================================================================
t_start <- Sys.time()
cat("[SF] loading readings ...\n")
d <- load_sf()
cat(sprintf("[SF] %s readings; generating synthetic sites ...\n",
            format(nrow(d$pts), big.mark = ",")))
sites <- gen_sites(d)
obs_hours <- sort(unique(hour(sites$mp$ts)))
saved <- .Random.seed
set.seed(HOUR_SEED_SF)
s_hour <- sample(obs_hours, N_SYN, replace = TRUE)
.Random.seed <- saved
sites$s_ts <- as.POSIXct(paste(sites$s_date), tz = "UTC") +
  s_hour * 3600 + sites$s_sec
cat(sprintf("[SF] %d observed collection hours (UTC): %s\n",
            length(obs_hours), paste(obs_hours, collapse = " ")))
cat("[SF] matching sweep ...\n")
sw <- run_sweep_idx(sites)
cat(sprintf("[SF] sweep done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

mp <- sites$mp
mp$lnx <- log(mp$nox)
cc0 <- which(CELLS$R == 2000 & CELLS$W == 72)   # baseline cell

# -- reproduction check against archived scripts/56 SF baseline ---------------
ib <- sw$idx_b[, cc0]; ic <- sw$idx_c[, cc0]
okb <- !is.na(ib); okc <- !is.na(ic); both <- okb & okc
cls_b <- match(mp$road_class[ib], cls_lab)
cm_nox <- mp %>% group_by(road_class) %>%
  summarise(m = mean(lnx), .groups = "drop")
mvec_nox <- setNames(cm_nox$m, cm_nox$road_class)
rep_chk <- c(
  match_b = round(100 * mean(okb), 1),
  match_c = round(100 * mean(okc), 1),
  fw      = round(100 * mean(cls_b[okb] == 4L), 1),
  own     = round(100 * mean(cls_b[okb] == match(sites$s_class, cls_lab)[okb]), 1),
  infl    = round(mean(mvec_nox[cls_lab[cls_b[okb]]]) -
                  mean(mvec_nox[sites$s_class[okb]]), 2),
  mn_b    = round(mean(mp$lnx[ib[both]]), 2),
  mn_c    = round(mean(mp$lnx[ic[both]]), 2),
  n_both  = sum(both))
arch_v <- c(ARCH$match_b, ARCH$match_c, ARCH$fw, ARCH$own, ARCH$infl,
            ARCH$mn_b, ARCH$mn_c, ARCH$n_both)
cat("reproduction check (this run vs archived scripts/56 SF baseline):\n")
print(rbind(this_run = rep_chk, archived = arch_v))
stopifnot(all(abs(rep_chk - arch_v) < 1e-8))
cat("REPRODUCED: scripts/56 SF baseline replicated bit-identically.\n")

# -- Part 1: SF PN support and class means ------------------------------------
pts_all <- d$pts                       # all SF readings (coords/ts valid)
n_all   <- nrow(pts_all)
pnv     <- pts_all$PN[is.finite(pts_all$PN)]
sh_miss <- 100 * (1 - length(pnv) / n_all)
sh_pn1  <- 100 * mean(is.finite(pts_all$PN1))
sh_neg  <- 100 * mean(pnv < 0)
sh_zero <- 100 * mean(pnv == 0)
sh_nonpos <- 100 * mean(pnv <= 0)
use_log <- sh_nonpos < 0.5
cls_pn_all <- pts_all %>% filter(is.finite(PN)) %>%
  group_by(road_class) %>%
  summarise(n = n(), mean_pn = mean(PN), med_pn = median(PN),
            mean_log = mean(log(PN[PN > 0])),
            sh_nonpos = 100 * mean(PN <= 0), .groups = "drop") %>%
  arrange(match(road_class, cls_lab))
# class means over the MATCHING POOL (valid NOx, as in scripts/56), used to
# value the composition gap so pool and valuation are internally consistent
cm_pn <- mp %>% filter(is.finite(PN)) %>% group_by(road_class) %>%
  summarise(m = mean(PN), .groups = "drop")
mvec_pn <- setNames(cm_pn$m, cm_pn$road_class)
# log-scale pool means (positive PN only; primary scale if use_log)
cm_lpn <- mp %>% filter(is.finite(PN), PN > 0) %>% group_by(road_class) %>%
  summarise(m = mean(log(PN)), .groups = "drop")
mvec_lpn <- setNames(cm_lpn$m, cm_lpn$road_class)

# -- Part 2: PN read off the identical matches --------------------------------
pn_cell <- function(cc) {
  ib <- sw$idx_b[, cc]; ic <- sw$idx_c[, cc]
  okb <- !is.na(ib)
  clsb <- match(mp$road_class[ib], cls_lab)
  infl <- mean(mvec_pn[cls_lab[clsb[okb]]]) - mean(mvec_pn[sites$s_class[okb]])
  infl_log <- mean(mvec_lpn[cls_lab[clsb[okb]]]) -
              mean(mvec_lpn[sites$s_class[okb]])
  both <- okb & !is.na(ic) & is.finite(mp$PN[ib]) & is.finite(mp$PN[ic])
  bothp <- both & mp$PN[ib] > 0 & mp$PN[ic] > 0
  tibble(label = CELLS$label[cc],
         infl_log = infl_log,
         gap_log = mean(log(mp$PN[ib[bothp]])) -
                   mean(log(mp$PN[ic[bothp]])),
         n_bothp = sum(bothp),
         infl_lvl = infl,
         mn_b = mean(mp$PN[ib[both]]),
         mn_c = mean(mp$PN[ic[both]]),
         gap_lvl = mean(mp$PN[ib[both]]) - mean(mp$PN[ic[both]]),
         n_both = sum(both),
         fw = 100 * mean(clsb[okb] == 4L))
}
pn_res <- map_dfr(seq_len(nrow(CELLS)), pn_cell)
b <- pn_res[cc0, ]
own_mean <- mean(mvec_pn[sites$s_class[!is.na(sw$idx_b[, cc0])]])
infl_pct <- 100 * b$infl_lvl / own_mean
gap_pct  <- 100 * b$gap_lvl / b$mn_c

# verdict (registered prediction: POSITIVE artifact, like NOx)
fw_gt_loc <- mvec_pn["freeway"] > mvec_pn["local"]
pos_infl  <- b$infl_log > 0
pos_gap   <- b$gap_log > 0
verdict <- if (fw_gt_loc && pos_infl && pos_gap) "CONFIRMED" else
  if (!fw_gt_loc && !pos_infl && !pos_gap) "NOT CONFIRMED" else "MIXED"

# -- write table --------------------------------------------------------------
out <- file("output/tables/ufp_replication.txt", open = "wt")
w <- function(...) cat(..., "\n", sep = "", file = out)
w("UFP (particle number) replication (registered prediction) -- scripts/65 -- ",
  format(Sys.Date()))
w("Prediction: ultrafine particles are a PRIMARY traffic pollutant, emitted")
w("directly by combustion, so PN should be HIGHEST on freeways and the")
w("siting-composition mechanism predicts a POSITIVE artifact -- blind joins")
w("OVERSTATE PN at non-freeway sites, same sign as NOx, mirror of O3")
w("(scripts/64).")
w("SF only: the CA release carries PN1-PN5 columns; the Oakland release does not.")
w("")
w("Design: the scripts/56 SF experiment replicated bit-identically (seed")
w("20260720; hour sub-stream 20260721; same 2,000 non-freeway sites, days,")
w("hours, matching pool = readings with valid NOx > 0, and combined-metric")
w("join). The JOIN IS UNCHANGED: each site receives exactly the reading it")
w("received in scripts/56; we read off that reading's PN instead of its NOx.")
w("Reproduction check vs archived scripts/56 SF baseline: PASSED (exact on")
w("all reported figures).")
w("")
w("== PN measure and support (all SF readings, n = ",
  format(n_all, big.mark = ","), ") ==")
w("  measure: TOTAL PN = PN1 + PN2 + PN3 + PN4 + PN5 (particle counts across")
w("  the release's five size bins; requires all five bins valid)")
w(sprintf("  coverage: total PN valid for %.1f%% of readings (PN1 alone %.1f%%,",
          100 - sh_miss, sh_pn1))
w("  so requiring all five bins costs essentially nothing; total PN chosen")
w("  as the natural particle-number concentration)")
w(sprintf("  among valid PN (n = %s): negative %.2f%%, zero %.2f%%, non-positive %.2f%%",
          format(length(pnv), big.mark = ","), sh_neg, sh_zero, sh_nonpos))
if (use_log) {
  w("  scale decision: non-positive share below 0.5%, so support permits")
  w("  logs; log(PN) reported as PRIMARY for comparability with the NOx")
  w("  log-point figures, levels (counts) alongside (non-positive readings")
  w("  dropped from log computations only).")
} else {
  w("  scale decision: LEVELS (counts) only -- non-positive share too large")
  w("  for logs without truncation.")
}
w("")
w("== SF mean PN by road class (all valid-PN readings) ==")
for (i in seq_len(nrow(cls_pn_all)))
  w(sprintf("  %-9s n = %9s  mean log(PN) %5.2f  mean %8s  median %8s  non-positive %5.2f%%",
            cls_pn_all$road_class[i],
            format(cls_pn_all$n[i], big.mark = ","),
            cls_pn_all$mean_log[i],
            format(round(cls_pn_all$mean_pn[i]), big.mark = ","),
            format(round(cls_pn_all$med_pn[i]), big.mark = ","),
            cls_pn_all$sh_nonpos[i]))
w(sprintf("  freeway > local? %s (pool means: %s; pool log means: %s)",
          ifelse(fw_gt_loc, "YES", "NO"),
          paste(sprintf("%s %s", names(mvec_pn),
                        format(round(mvec_pn), big.mark = ",")), collapse = ", "),
          paste(sprintf("%s %.2f", names(mvec_lpn), mvec_lpn), collapse = ", ")))
w("")
w("== Synthetic experiment, PN read off the identical scripts/56 matches ==")
w("   (baseline cell 2 km / +/-72 h; composition gap valued at matching-pool")
w("    class means; paired gap on the both-matched subset with valid PN at")
w("    both matched readings)")
w(sprintf("  blind freeway-assigned share (unchanged by construction): %.1f%%",
          b$fw))
w(sprintf("  composition-implied gap: %+.2f log points (levels: %+s counts, %+.1f%% of own-class mean %s)",
          b$infl_log, format(round(b$infl_lvl), big.mark = ","), infl_pct,
          format(round(own_mean), big.mark = ",")))
w(sprintf("  paired blind-vs-constrained gap (n = %d): %+.2f log points (levels: blind %s, constrained %s, gap %+s counts, %+.1f%%)",
          b$n_bothp, b$gap_log,
          format(round(b$mn_b), big.mark = ","),
          format(round(b$mn_c), big.mark = ","),
          format(round(b$gap_lvl), big.mark = ","), gap_pct))
w("  (n reflects PN instrument coverage: both matched readings must carry")
w("   valid PN in all five bins; the joins themselves are the scripts/56 joins.)")
w("")
w("  parameter sweep (same cells as scripts/55/56):")
w(sprintf("  %-18s %14s %15s %10s %10s %18s", "cell", "compo-infl log",
          "paired-log (n)", "blind PN", "constr PN", "paired-lvl (n)"))
for (i in seq_len(nrow(pn_res)))
  w(sprintf("  %-18s %+13.2f %+9.2f (%4d) %10s %10s %+11s (%d)",
            pn_res$label[i], pn_res$infl_log[i], pn_res$gap_log[i],
            pn_res$n_bothp[i],
            format(round(pn_res$mn_b[i]), big.mark = ","),
            format(round(pn_res$mn_c[i]), big.mark = ","),
            format(round(pn_res$gap_lvl[i]), big.mark = ","),
            pn_res$n_both[i]))
w("")
w("== VERDICT ==")
w(sprintf("  freeway mean above local mean: %s", ifelse(fw_gt_loc, "yes", "no")))
w(sprintf("  composition-implied gap positive: %s (%+.2f log points)",
          ifelse(pos_infl, "yes", "no"), b$infl_log))
w(sprintf("  paired blind-vs-constrained gap positive: %s (%+.2f log points)",
          ifelse(pos_gap, "yes", "no"), b$gap_log))
w(sprintf("  Registered prediction (positive, NOx-signed artifact): %s", verdict))
w("  (Comparison, same sites/joins: NOx composition inflation +0.81 log")
w("   points, paired gap +0.15 log points -- scripts/56; O3 composition gap")
w("   -0.18 log points [-3.96 ppb] -- scripts/64.)")
close(out)
cat(readLines("output/tables/ufp_replication.txt"), sep = "\n")
