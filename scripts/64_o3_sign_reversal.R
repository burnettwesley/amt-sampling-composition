# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 64_o3_sign_reversal.R
# =====================
# Registered-prediction test: the siting-composition mechanism predicts a
# SIGN-REVERSED artifact for ozone. Fresh NO titrates O3 near roads
# (NO + O3 -> NO2 + O2), so O3 is LOW where NOx is high; a blind join that
# imports major-road readings should therefore UNDERSTATE O3 at non-freeway
# sites (negative composition inflation), the mirror image of the NOx
# overstatement documented by scripts/45/55/56.
#
# DESIGN: SF only (the California_201605_201709 release carries an O3 column;
# the Oakland release does not). We replicate the scripts/56 SF synthetic
# experiment BIT-IDENTICALLY -- same seed (20260720, SF processed first, so
# the RNG stream needs no Oakland pass), same 2,000 non-freeway sites, same
# collection days, same uniform hour draws from the isolated sub-stream
# (seed 20260721), same discarded burn draw, same matching pool (readings
# with valid NOx > 0) and the same combined-metric join. The JOIN IS
# UNCHANGED: each synthetic site is matched to exactly the reading it
# received in scripts/56; we then read off that reading's O3 instead of its
# NOx. A hard reproduction check confirms the baseline-cell NOx numbers
# equal the archived scripts/56 values before any O3 result is written.
#
# O3 SUPPORT (handled honestly, see output): O3 is reported in ppb. Missing
# O3 is common (57% of SF readings; the O3 instrument was not always
# aboard), but among valid readings non-positive values are essentially
# absent (0.00%; a handful of freeway zeros), so the support permits logs.
# We therefore report LEVELS (ppb) as the primary scale and log(O3)
# alongside, for direct comparability with the NOx log-point figures; if
# the non-positive share had exceeded 0.5% the script would fall back to
# levels only, since dropping titrated near-zero readings would truncate
# exactly the signal under test. The script documents the actual shares.
#
# NO outcome data of any kind are read (no crashes, no violations).
#
# OUTPUT: output/tables/o3_sign_reversal.txt
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

# -- loading (scripts/56 verbatim, plus the O3 column) ------------------------
load_sf <- function() {
  pts <- read_csv(CA_FILE,
                  col_select = c(Date_Time, Latitude, Longitude, O3, NO2, NO),
                  show_col_types = FALSE) %>%
    mutate(ts = ymd_hms(Date_Time, tz = "UTC"),
           O3  = suppressWarnings(as.numeric(O3)),
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
  xy <- pts %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(UTM10N) %>% st_coordinates()
  pts$x <- xy[, 1]; pts$y <- xy[, 2]
  pts <- pts %>% mutate(nox = NO2 + NO, date = as.Date(ts))
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
       mp = pts %>% select(ts, x, y, road_class, nox, O3))
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

# -- Part 1: SF O3 support and class means ------------------------------------
pts_all <- d$pts                       # all SF readings (coords/ts valid)
n_all   <- nrow(pts_all)
o3v     <- pts_all$O3[is.finite(pts_all$O3)]
sh_miss <- 100 * (1 - length(o3v) / n_all)
sh_neg  <- 100 * mean(o3v < 0)
sh_zero <- 100 * mean(o3v == 0)
sh_nonpos <- 100 * mean(o3v <= 0)
use_log <- sh_nonpos < 0.5
cls_o3_all <- pts_all %>% filter(is.finite(O3)) %>%
  group_by(road_class) %>%
  summarise(n = n(), mean_ppb = mean(O3), med_ppb = median(O3),
            sh_nonpos = 100 * mean(O3 <= 0), .groups = "drop") %>%
  arrange(match(road_class, cls_lab))
# class means over the MATCHING POOL (valid NOx, as in scripts/56), used to
# value the composition gap so pool and valuation are internally consistent
cm_o3 <- mp %>% filter(is.finite(O3)) %>% group_by(road_class) %>%
  summarise(m = mean(O3), .groups = "drop")
mvec_o3 <- setNames(cm_o3$m, cm_o3$road_class)
# log-scale pool means (positive O3 only; used only if use_log)
cm_lo3 <- mp %>% filter(is.finite(O3), O3 > 0) %>% group_by(road_class) %>%
  summarise(m = mean(log(O3)), .groups = "drop")
mvec_lo3 <- setNames(cm_lo3$m, cm_lo3$road_class)

# -- Part 2: O3 read off the identical matches --------------------------------
o3_cell <- function(cc) {
  ib <- sw$idx_b[, cc]; ic <- sw$idx_c[, cc]
  okb <- !is.na(ib)
  clsb <- match(mp$road_class[ib], cls_lab)
  infl <- mean(mvec_o3[cls_lab[clsb[okb]]]) - mean(mvec_o3[sites$s_class[okb]])
  both <- okb & !is.na(ic) & is.finite(mp$O3[ib]) & is.finite(mp$O3[ic])
  infl_log <- mean(mvec_lo3[cls_lab[clsb[okb]]]) -
              mean(mvec_lo3[sites$s_class[okb]])
  bothp <- both & mp$O3[ib] > 0 & mp$O3[ic] > 0
  tibble(label = CELLS$label[cc],
         infl_ppb = infl,
         mn_b = mean(mp$O3[ib[both]]),
         mn_c = mean(mp$O3[ic[both]]),
         gap_ppb = mean(mp$O3[ib[both]]) - mean(mp$O3[ic[both]]),
         n_both = sum(both),
         infl_log = infl_log,
         gap_log = mean(log(mp$O3[ib[bothp]])) -
                   mean(log(mp$O3[ic[bothp]])),
         n_bothp = sum(bothp),
         fw = 100 * mean(clsb[okb] == 4L))
}
o3_res <- map_dfr(seq_len(nrow(CELLS)), o3_cell)
b <- o3_res[cc0, ]
own_mean <- mean(mvec_o3[sites$s_class[!is.na(sw$idx_b[, cc0])]])
infl_pct <- 100 * b$infl_ppb / own_mean
gap_pct  <- 100 * b$gap_ppb / b$mn_c

# verdict
fw_lt_loc <- mvec_o3["freeway"] < mvec_o3["local"]
neg_infl  <- b$infl_ppb < 0
neg_gap   <- b$gap_ppb < 0
verdict <- if (fw_lt_loc && neg_infl && neg_gap) "CONFIRMED" else
  if (!fw_lt_loc && !neg_infl && !neg_gap) "NOT CONFIRMED" else "MIXED"

# -- write table --------------------------------------------------------------
out <- file("output/tables/o3_sign_reversal.txt", open = "wt")
w <- function(...) cat(..., "\n", sep = "", file = out)
w("O3 sign-reversal test (registered prediction) -- scripts/64 -- ",
  format(Sys.Date()))
w("Prediction: fresh NO titrates O3 near roads (NO + O3 -> NO2 + O2), so O3")
w("is LOW where NOx is high; the siting-composition mechanism therefore")
w("predicts that blind joins UNDERSTATE O3 at non-freeway sites (negative")
w("composition inflation), the mirror of the NOx overstatement.")
w("SF only: the CA release carries an O3 column; the Oakland release does not.")
w("")
w("Design: the scripts/56 SF experiment replicated bit-identically (seed")
w("20260720; hour sub-stream 20260721; same 2,000 non-freeway sites, days,")
w("hours, matching pool = readings with valid NOx > 0, and combined-metric")
w("join). The JOIN IS UNCHANGED: each site receives exactly the reading it")
w("received in scripts/56; we read off that reading's O3 instead of its NOx.")
w("Reproduction check vs archived scripts/56 SF baseline: PASSED (exact on")
w("all reported figures).")
w("")
w("== O3 support (all SF readings, n = ", format(n_all, big.mark = ","),
  ") ==")
w(sprintf("  units: ppb (per the Aclima data release)"))
w(sprintf("  O3 missing: %.1f%% of readings (instrument not always aboard)",
          sh_miss))
w(sprintf("  among valid O3 (n = %s): negative %.2f%%, zero %.2f%%, non-positive %.2f%%",
          format(length(o3v), big.mark = ","), sh_neg, sh_zero, sh_nonpos))
if (use_log) {
  w("  scale decision: non-positive share below 0.5%, so support permits")
  w("  logs; levels (ppb) reported as primary, log(O3) alongside for")
  w("  comparability with the NOx log-point figures (non-positive readings")
  w("  dropped from log computations only).")
} else {
  w("  scale decision: LEVELS (ppb) only -- log(O3) not used, because")
  w("  dropping non-positive O3 would truncate exactly the titrated")
  w("  near-road readings whose depletion the prediction concerns.")
}
w("")
w("== SF mean O3 by road class (all valid-O3 readings) ==")
for (i in seq_len(nrow(cls_o3_all)))
  w(sprintf("  %-9s n = %9s  mean %6.2f ppb  median %6.2f ppb  non-positive %5.2f%%",
            cls_o3_all$road_class[i],
            format(cls_o3_all$n[i], big.mark = ","),
            cls_o3_all$mean_ppb[i], cls_o3_all$med_ppb[i],
            cls_o3_all$sh_nonpos[i]))
w(sprintf("  freeway < local? %s (freeway %.2f vs local %.2f ppb; matching-pool means: %s)",
          ifelse(fw_lt_loc, "YES", "NO"),
          cls_o3_all$mean_ppb[cls_o3_all$road_class == "freeway"],
          cls_o3_all$mean_ppb[cls_o3_all$road_class == "local"],
          paste(sprintf("%s %.2f", names(mvec_o3), mvec_o3), collapse = ", ")))
w("")
w("== Synthetic experiment, O3 read off the identical scripts/56 matches ==")
w("   (baseline cell 2 km / +/-72 h; composition gap valued at matching-pool")
w("    class means of O3 in ppb; paired gap on the both-matched subset with")
w("    finite O3 at both matched readings)")
w(sprintf("  blind freeway-assigned share (unchanged by construction): %.1f%%",
          b$fw))
w(sprintf("  composition-implied gap: %+.2f ppb (%+.1f%% of the own-class mean %.2f ppb)",
          b$infl_ppb, infl_pct, own_mean))
w(sprintf("  paired blind-vs-constrained gap (n = %d): blind %.2f, constrained %.2f (gap %+.2f ppb, %+.1f%%)",
          b$n_both, b$mn_b, b$mn_c, b$gap_ppb, gap_pct))
w("  (n reflects O3 instrument coverage: both matched readings must carry")
w("   a valid O3 value; the joins themselves are the scripts/56 joins.)")
if (use_log) {
  w(sprintf("  log scale: composition-implied gap %+.2f log points; paired gap %+.2f log points (n = %d)",
            b$infl_log, b$gap_log, b$n_bothp))
  w(sprintf("  pool class means of log(O3): %s",
            paste(sprintf("%s %.2f", names(mvec_lo3), mvec_lo3),
                  collapse = ", ")))
}
w("")
w("  parameter sweep (same cells as scripts/55/56):")
w(sprintf("  %-18s %14s %10s %10s %16s", "cell", "compo-infl ppb",
          "blind O3", "constr O3", "paired-gap (n)"))
for (i in seq_len(nrow(o3_res)))
  w(sprintf("  %-18s %+13.2f %10.2f %10.2f %+9.2f (%d)",
            o3_res$label[i], o3_res$infl_ppb[i], o3_res$mn_b[i],
            o3_res$mn_c[i], o3_res$gap_ppb[i], o3_res$n_both[i]))
w("")
w("== VERDICT ==")
w(sprintf("  freeway mean below local mean: %s", ifelse(fw_lt_loc, "yes", "no")))
w(sprintf("  composition-implied gap negative: %s (%+.2f ppb)",
          ifelse(neg_infl, "yes", "no"), b$infl_ppb))
w(sprintf("  paired blind-vs-constrained gap negative: %s (%+.2f ppb)",
          ifelse(neg_gap, "yes", "no"), b$gap_ppb))
w(sprintf("  Sign-reversal prediction: %s", verdict))
w("  (NOx comparison, same sites/joins: composition inflation +0.81 log")
w("   points, paired gap +0.15 log points -- scripts/56.)")
close(out)
cat(readLines("output/tables/o3_sign_reversal.txt"), sep = "\n")
