# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 67_amt_bootstrap.R
# ==================
# Uncertainty quantification for the synthetic-outcome experiment: a
# nonparametric bootstrap over synthetic SITES, the experiment's sampling
# unit. B = 999 resamples of the 2,000 sites per city (with replacement),
# recomputing per resample:
#   (a) the composition-implied inflation -- both cities for NOx; for SF also
#       O3 (log scale, archived -0.18) and UFP/PN (log scale, archived +0.11),
#       read off the IDENTICAL scripts/56 matches per scripts/64/65;
#   (b) the paired blind-vs-constrained gap for each pollutant.
# Percentile 95% CIs (quantile type 7, probs .025/.975). Because all three SF
# pollutants are read off the same matches, the SAME resample indices are used
# across pollutants, so pairwise composition-gap differences (NOx-UFP,
# UFP-O3) are paired within resample and their CIs answer whether the three
# artifacts are mutually distinct.
#
# Also: a bootstrap CI for the scripts/42 cross-car correlation, resampling
# the 1,198 dual-coverage cell-hours (recomputed here from the Oakland raw
# file, which this script loads anyway; the scripts/42 panel construction is
# replicated and the n = 1,198 count is asserted before any CI is written).
#
# Class-mean valuation vectors (matching-pool class means of log NOx, O3 ppb,
# log O3, PN, log PN) are held FIXED at their full-pool values: they average
# hundreds of thousands to millions of readings, so their sampling error is
# negligible next to the site resampling, and fixing them keeps each
# resample's inflation valued on the same scale as the archived point
# estimates.
#
# RNG DISCIPLINE: site construction replicates scripts/56 bit-identically
# (seed 20260720, SF first; hour sub-streams 20260721/20260722; discarded
# burn draw). Hard reproduction checks against the archived scripts/56/64/65
# baselines run BEFORE any bootstrap. The bootstrap itself uses a DIFFERENT,
# documented seed (20260805), set only AFTER all site construction and
# matching are complete, so the archived baselines are untouched.
#
# NO outcome data of any kind are read (no crashes, no violations).
#
# OUTPUT: output/tables/amt_bootstrap.txt
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(lubridate)
})
if (!file.exists("scripts/45_sf_siting.R"))
  stop("Run this script from the repository root.")

set.seed(20260720)   # same seed as scripts/45/56; SF is processed first
HOUR_SEED <- c(SF = 20260721, Oakland = 20260722)
BOOT_SEED <- 20260805   # bootstrap-only seed, set AFTER site construction
B         <- 999

UTM10N  <- 32610
N_SYN   <- 2000
N_BURN  <- 1217L     # length of scripts/45's crash-hour vector (burn only)
CA_FILE <- Sys.getenv("AQ_CA_FILE", "data/raw/aclima/California_201605_201709_GoogleAclimaAQ.txt")
SFBB    <- c(lat_min = 37.708, lat_max = 37.833, lon_min = -122.515, lon_max = -122.355)
cls_lab <- c("local", "collector", "arterial", "freeway")

# baseline cell only (2 km / +/-72 h); restricting the candidate window to
# the baseline caps leaves the baseline-cell matches identical to scripts/56
R_CAP <- 2000; W_CAP <- 72

# archived baselines (scripts/56 both cities; scripts/64 O3; scripts/65 PN)
ARCH56 <- list(
  SF      = c(match_b = 41.9, match_c = 34.3, fw = 36.0, own = 21.6,
              infl = 0.81, mn_b = 2.71, mn_c = 2.56, n_both = 686),
  Oakland = c(match_b = 44.7, match_c = 30.1, fw = 64.7, own = 11.6,
              infl = 1.05, mn_b = 3.74, mn_c = 3.45, n_both = 602))
ARCH64 <- c(infl_ppb = -3.96, gap_ppb = -0.40, n_both = 182,
            infl_log = -0.18, gap_log = 0.00)
ARCH65 <- c(infl_log = 0.11, gap_log = 0.01, n_bothp = 615)
ARCH42 <- c(n_cc = 1198, corr = 0.656)

# -- loading (scripts/56 verbatim; SF adds O3 + PN cols, Oakland adds car id) --
load_city <- function(which) {
  if (which == "SF") {
    pts <- read_csv(CA_FILE,
                    col_select = c(Date_Time, Latitude, Longitude,
                                   O3, PN1, PN2, PN3, PN4, PN5, NO2, NO),
                    show_col_types = FALSE) %>%
      mutate(ts = ymd_hms(Date_Time, tz = "UTC"),
             across(c(O3, PN1, PN2, PN3, PN4, PN5, NO2, NO),
                    ~ suppressWarnings(as.numeric(.x)))) %>%
      filter(!is.na(Latitude), !is.na(Longitude), !is.na(ts),
             Latitude >= SFBB["lat_min"], Latitude <= SFBB["lat_max"],
             Longitude >= SFBB["lon_min"], Longitude <= SFBB["lon_max"])
    rc <- read_csv("data/processed/sf_reading_road_class.csv.gz",
                   show_col_types = FALSE)
    stopifnot(nrow(rc) == nrow(pts))
    pts$road_class <- rc$road_class
    pts <- pts %>% mutate(PN = PN1 + PN2 + PN3 + PN4 + PN5)  # NA if any bin NA
    roads <- st_read("data/processed/sf_osm_roads.gpkg", quiet = TRUE)
  } else {
    pts <- read_csv("data/raw/aclima/Oakland_201505_201605_GoogleAclimaAQ.txt",
                    col_select = c(Date_Time, Car_Identifier, Latitude,
                                   Longitude, NO2, NO),
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
       mp = pts %>% select(any_of(c("ts", "x", "y", "road_class",
                                    "nox", "O3", "PN"))))
}

# -- baseline-cell join (scripts/64/65 verbatim, single cell), returns indices -
match_baseline <- function(sites) {
  mp <- sites$mp
  orig <- seq_len(nrow(mp))
  o <- order(as.numeric(mp$ts))
  tsn <- as.numeric(mp$ts)[o]
  px <- mp$x[o]; py <- mp$y[o]
  cint <- match(mp$road_class, cls_lab)[o]
  orig <- orig[o]
  s_cint <- match(sites$s_class, cls_lab)
  t0v <- as.numeric(sites$s_ts)
  idx_b <- rep(NA_integer_, N_SYN); idx_c <- idx_b
  for (i in seq_len(N_SYN)) {
    t0 <- t0v[i]
    lo <- findInterval(t0 - W_CAP * 3600, tsn, left.open = TRUE) + 1L
    hi <- findInterval(t0 + W_CAP * 3600, tsn)
    if (hi < lo) next
    ii <- lo:hi
    dx <- px[ii] - sites$sxy[i, 1]; dy <- py[ii] - sites$sxy[i, 2]
    sp <- sqrt(dx * dx + dy * dy)
    k <- sp <= R_CAP
    if (!any(k)) next
    spk <- sp[k]
    tk  <- abs(tsn[ii][k] - t0) / 3600
    ck  <- cint[ii][k]; ok <- orig[ii][k]
    sel <- which(spk <= R_CAP & tk <= W_CAP)
    if (!length(sel)) next
    dd <- sqrt((spk[sel] / 500)^2 + (tk[sel] / 6)^2)
    tie <- sel[dd == min(dd)]
    j <- tie[which.min(ok[tie])]
    idx_b[i] <- ok[j]
    selc <- sel[ck[sel] == s_cint[i]]
    if (length(selc)) {
      ddc <- sqrt((spk[selc] / 500)^2 + (tk[selc] / 6)^2)
      tiec <- selc[ddc == min(ddc)]
      jc <- tiec[which.min(ok[tiec])]
      idx_c[i] <- ok[jc]
    }
  }
  list(idx_b = idx_b, idx_c = idx_c)
}

# =============================================================================
# Pass 1: construct sites and matches for both cities (RNG stream = scripts/56)
# =============================================================================
cache <- list(); cc_panel <- NULL
for (city in c("SF", "Oakland")) {
  t_start <- Sys.time()
  cat(sprintf("[%s] loading readings ...\n", city))
  d <- load_city(city)
  cat(sprintf("[%s] %s readings; generating synthetic sites ...\n",
              city, format(nrow(d$pts), big.mark = ",")))
  sites <- gen_sites(d)
  obs_hours <- sort(unique(hour(sites$mp$ts)))
  saved <- .Random.seed
  set.seed(HOUR_SEED[[city]])
  s_hour <- sample(obs_hours, N_SYN, replace = TRUE)
  .Random.seed <- saved
  sites$s_ts <- as.POSIXct(paste(sites$s_date), tz = "UTC") +
    s_hour * 3600 + sites$s_sec
  cat(sprintf("[%s] baseline-cell join ...\n", city))
  mm <- match_baseline(sites)
  mp <- sites$mp
  mp$lnx <- log(mp$nox)

  # pool class-mean valuation vectors (fixed; see header)
  cmv <- function(v) {
    tb <- tibble(road_class = mp$road_class, v = v) %>%
      filter(is.finite(v)) %>% group_by(road_class) %>%
      summarise(m = mean(v), .groups = "drop")
    setNames(tb$m, tb$road_class)
  }
  mvec_nox <- cmv(mp$lnx)

  ib <- mm$idx_b; ic <- mm$idx_c
  okb <- !is.na(ib); okc <- !is.na(ic); both <- okb & okc
  cls_b <- match(mp$road_class[ib], cls_lab)

  # reproduction check vs archived scripts/56 baseline
  chk <- c(match_b = round(100 * mean(okb), 1),
           match_c = round(100 * mean(okc), 1),
           fw      = round(100 * mean(cls_b[okb] == 4L), 1),
           own     = round(100 * mean(cls_b[okb] ==
                             match(sites$s_class, cls_lab)[okb]), 1),
           infl    = round(mean(mvec_nox[cls_lab[cls_b[okb]]]) -
                           mean(mvec_nox[sites$s_class[okb]]), 2),
           mn_b    = round(mean(mp$lnx[ib[both]]), 2),
           mn_c    = round(mean(mp$lnx[ic[both]]), 2),
           n_both  = sum(both))
  cat(sprintf("[%s] reproduction check vs archived scripts/56:\n", city))
  print(rbind(this_run = chk, archived = ARCH56[[city]]))
  stopifnot(all(abs(chk - ARCH56[[city]]) < 1e-8))
  cat(sprintf("[%s] REPRODUCED scripts/56 baseline.\n", city))

  # per-site cached quantities (NA outside the relevant subset)
  d_nox <- ifelse(okb, mvec_nox[cls_lab[cls_b]] -
                       mvec_nox[sites$s_class], NA_real_)
  g_nox <- ifelse(both, mp$lnx[ib] - mp$lnx[ic], NA_real_)
  cty <- list(okb = okb, d_nox = d_nox, g_nox = g_nox)

  if (city == "SF") {
    mvec_o3  <- cmv(mp$O3)
    mvec_lo3 <- cmv(ifelse(is.finite(mp$O3) & mp$O3 > 0, log(mp$O3), NA_real_))
    mvec_lpn <- cmv(ifelse(is.finite(mp$PN) & mp$PN > 0, log(mp$PN), NA_real_))
    both_o3  <- both & is.finite(mp$O3[ib]) & is.finite(mp$O3[ic])
    bothp_o3 <- both_o3 & mp$O3[ib] > 0 & mp$O3[ic] > 0
    both_pn  <- both & is.finite(mp$PN[ib]) & is.finite(mp$PN[ic])
    bothp_pn <- both_pn & mp$PN[ib] > 0 & mp$PN[ic] > 0
    cty$d_o3_log  <- ifelse(okb, mvec_lo3[cls_lab[cls_b]] -
                                 mvec_lo3[sites$s_class], NA_real_)
    cty$d_o3_ppb  <- ifelse(okb, mvec_o3[cls_lab[cls_b]] -
                                 mvec_o3[sites$s_class], NA_real_)
    cty$g_o3_ppb  <- ifelse(both_o3, mp$O3[ib] - mp$O3[ic], NA_real_)
    cty$g_o3_log  <- ifelse(bothp_o3, log(mp$O3[ib]) - log(mp$O3[ic]), NA_real_)
    cty$d_pn_log  <- ifelse(okb, mvec_lpn[cls_lab[cls_b]] -
                                 mvec_lpn[sites$s_class], NA_real_)
    cty$g_pn_log  <- ifelse(bothp_pn, log(mp$PN[ib]) - log(mp$PN[ic]), NA_real_)
    # reproduction checks vs archived scripts/64 and /65 baselines
    chk64 <- c(infl_ppb = round(mean(cty$d_o3_ppb, na.rm = TRUE), 2),
               gap_ppb  = round(mean(cty$g_o3_ppb, na.rm = TRUE), 2),
               n_both   = sum(both_o3),
               infl_log = round(mean(cty$d_o3_log, na.rm = TRUE), 2),
               gap_log  = round(mean(cty$g_o3_log, na.rm = TRUE), 2))
    cat("[SF] reproduction check vs archived scripts/64 (O3):\n")
    print(rbind(this_run = chk64, archived = ARCH64))
    stopifnot(all(abs(chk64 - ARCH64) < 1e-8))
    chk65 <- c(infl_log = round(mean(cty$d_pn_log, na.rm = TRUE), 2),
               gap_log  = round(mean(cty$g_pn_log, na.rm = TRUE), 2),
               n_bothp  = sum(bothp_pn))
    cat("[SF] reproduction check vs archived scripts/65 (PN):\n")
    print(rbind(this_run = chk65, archived = ARCH65))
    stopifnot(all(abs(chk65 - ARCH65) < 1e-8))
    cat("[SF] REPRODUCED scripts/64 and /65 baselines.\n")
  } else {
    # cross-car dual-coverage panel, scripts/42 construction verbatim
    cc_panel <- d$pts %>%
      mutate(cell = paste(floor(Longitude / 0.0012),
                          floor(Latitude / 0.0009)),
             hr   = floor_date(ts, "hour")) %>%
      filter(!is.na(nox), nox > 0) %>%
      group_by(cell, hr, Car_Identifier) %>%
      summarise(n = n(), mean_nox = mean(nox), .groups = "drop") %>%
      filter(n >= 5) %>%
      pivot_wider(id_cols = c(cell, hr), names_from = Car_Identifier,
                  values_from = mean_nox) %>%
      filter(!is.na(Car_A), !is.na(Car_B)) %>%
      mutate(lA = log(Car_A), lB = log(Car_B))
    cat(sprintf("[Oakland] cross-car panel: %d dual-coverage cell-hours, corr %.3f\n",
                nrow(cc_panel), cor(cc_panel$lA, cc_panel$lB)))
    stopifnot(nrow(cc_panel) == ARCH42["n_cc"],
              abs(round(cor(cc_panel$lA, cc_panel$lB), 3) - ARCH42["corr"]) < 1e-8)
    cat("[Oakland] REPRODUCED scripts/42 panel (n = 1,198, corr 0.656).\n")
  }
  cache[[city]] <- cty
  cat(sprintf("[%s] done in %.1f min\n", city,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  rm(d, sites, mp); invisible(gc())
}

# =============================================================================
# Pass 2: bootstrap (fresh, documented seed; archived baselines untouched)
# =============================================================================
set.seed(BOOT_SEED)
boot_stat <- function(v, idx) {           # mean of v over resampled sites,
  colMeans(matrix(v[idx], nrow(idx), B),  # NA-dropping = subset ratio est.
           na.rm = TRUE)
}
ci <- function(x) quantile(x, c(.025, .975), names = FALSE)  # type 7
fmt_ci <- function(pt, x, dg = 2)
  sprintf(paste0("%+.", dg, "f [%+.", dg, "f, %+.", dg, "f]"),
          pt, ci(x)[1], ci(x)[2])
excl0 <- function(x) { q <- ci(x); ifelse(q[1] > 0 | q[2] < 0, "YES", "NO") }

idx_SF  <- matrix(sample.int(N_SYN, N_SYN * B, replace = TRUE), N_SYN, B)
idx_OK  <- matrix(sample.int(N_SYN, N_SYN * B, replace = TRUE), N_SYN, B)
idx_CC  <- matrix(sample.int(nrow(cc_panel), nrow(cc_panel) * B,
                             replace = TRUE), nrow(cc_panel), B)

sf <- cache$SF; ok <- cache$Oakland
bs <- list(
  infl_nox_sf = boot_stat(sf$d_nox,    idx_SF),
  gap_nox_sf  = boot_stat(sf$g_nox,    idx_SF),
  infl_nox_ok = boot_stat(ok$d_nox,    idx_OK),
  gap_nox_ok  = boot_stat(ok$g_nox,    idx_OK),
  infl_o3     = boot_stat(sf$d_o3_log, idx_SF),
  infl_o3_ppb = boot_stat(sf$d_o3_ppb, idx_SF),
  gap_o3_ppb  = boot_stat(sf$g_o3_ppb, idx_SF),
  gap_o3_log  = boot_stat(sf$g_o3_log, idx_SF),
  infl_pn     = boot_stat(sf$d_pn_log, idx_SF),
  gap_pn_log  = boot_stat(sf$g_pn_log, idx_SF))
# paired pollutant differences (same resample indices within SF)
bs$diff_nox_pn <- bs$infl_nox_sf - bs$infl_pn
bs$diff_pn_o3  <- bs$infl_pn    - bs$infl_o3
bs$diff_nox_o3 <- bs$infl_nox_sf - bs$infl_o3
# cross-car correlation
bs$cc_corr <- vapply(seq_len(B), function(b) {
  i <- idx_CC[, b]; cor(cc_panel$lA[i], cc_panel$lB[i])
}, numeric(1))

pt <- list(
  infl_nox_sf = mean(sf$d_nox, na.rm = TRUE),
  gap_nox_sf  = mean(sf$g_nox, na.rm = TRUE),
  infl_nox_ok = mean(ok$d_nox, na.rm = TRUE),
  gap_nox_ok  = mean(ok$g_nox, na.rm = TRUE),
  infl_o3     = mean(sf$d_o3_log, na.rm = TRUE),
  infl_o3_ppb = mean(sf$d_o3_ppb, na.rm = TRUE),
  gap_o3_ppb  = mean(sf$g_o3_ppb, na.rm = TRUE),
  gap_o3_log  = mean(sf$g_o3_log, na.rm = TRUE),
  infl_pn     = mean(sf$d_pn_log, na.rm = TRUE),
  gap_pn_log  = mean(sf$g_pn_log, na.rm = TRUE),
  cc_corr     = cor(cc_panel$lA, cc_panel$lB))
pt$diff_nox_pn <- pt$infl_nox_sf - pt$infl_pn
pt$diff_pn_o3  <- pt$infl_pn    - pt$infl_o3
pt$diff_nox_o3 <- pt$infl_nox_sf - pt$infl_o3

# =============================================================================
# write table
# =============================================================================
out <- file("output/tables/amt_bootstrap.txt", open = "wt")
w <- function(...) cat(..., "\n", sep = "", file = out)
w("AMT uncertainty quantification: site bootstrap -- scripts/67 -- ",
  format(Sys.Date()))
w("Nonparametric bootstrap over synthetic SITES (the sampling unit):")
w(sprintf("B = %d resamples of the %d sites per city, with replacement,", B, N_SYN))
w("recomputing the composition-implied inflation and the paired")
w("blind-vs-constrained gap per resample; percentile 95% CIs (probs")
w(".025/.975, quantile type 7). Site construction and joins replicate")
w("scripts/56/64/65 bit-identically (seed 20260720; hour sub-streams")
w("20260721/20260722); hard reproduction checks PASSED for the scripts/56")
w("SF and Oakland baselines, the scripts/64 O3 figures, the scripts/65 PN")
w(sprintf("figures, and the scripts/42 panel before any CI was computed."))
w(sprintf("Bootstrap seed: %d, set AFTER all site construction, so archived", BOOT_SEED))
w("baselines are untouched. Class-mean valuation vectors held fixed at")
w("full-pool values (sampling error negligible at pool n; see script header).")
w("The three SF pollutants are read off the identical matches and share")
w("resample indices, so pairwise differences are paired within resample.")
w("")
w("== Composition-implied inflation (blind class mix vs own-class mix) ==")
w(sprintf("  Oakland NOx (log points): %s", fmt_ci(pt$infl_nox_ok, bs$infl_nox_ok)))
w(sprintf("  SF NOx      (log points): %s", fmt_ci(pt$infl_nox_sf, bs$infl_nox_sf)))
w(sprintf("  SF O3       (log points): %s  excludes zero: %s",
          fmt_ci(pt$infl_o3, bs$infl_o3), excl0(bs$infl_o3)))
w(sprintf("  SF O3       (ppb):        %s", fmt_ci(pt$infl_o3_ppb, bs$infl_o3_ppb)))
w(sprintf("  SF UFP/PN   (log points): %s  excludes zero: %s",
          fmt_ci(pt$infl_pn, bs$infl_pn), excl0(bs$infl_pn)))
w("")
w("== Paired blind-vs-constrained gap (both-matched subset) ==")
w(sprintf("  Oakland NOx (log points): %s", fmt_ci(pt$gap_nox_ok, bs$gap_nox_ok)))
w(sprintf("  SF NOx      (log points): %s", fmt_ci(pt$gap_nox_sf, bs$gap_nox_sf)))
w(sprintf("  SF O3       (ppb):        %s  excludes zero: %s",
          fmt_ci(pt$gap_o3_ppb, bs$gap_o3_ppb), excl0(bs$gap_o3_ppb)))
w(sprintf("  SF O3       (log points): %s  excludes zero: %s",
          fmt_ci(pt$gap_o3_log, bs$gap_o3_log), excl0(bs$gap_o3_log)))
w(sprintf("  SF UFP/PN   (log points): %s  excludes zero: %s",
          fmt_ci(pt$gap_pn_log, bs$gap_pn_log), excl0(bs$gap_pn_log)))
w("")
w("== Pairwise composition-gap differences (log points, SF; paired ==")
w("   within resample -- identical sites and joins across pollutants)")
w(sprintf("  NOx - UFP: %s  excludes zero: %s",
          fmt_ci(pt$diff_nox_pn, bs$diff_nox_pn), excl0(bs$diff_nox_pn)))
w(sprintf("  UFP - O3:  %s  excludes zero: %s",
          fmt_ci(pt$diff_pn_o3, bs$diff_pn_o3), excl0(bs$diff_pn_o3)))
w(sprintf("  NOx - O3:  %s  excludes zero: %s",
          fmt_ci(pt$diff_nox_o3, bs$diff_nox_o3), excl0(bs$diff_nox_o3)))
w("")
w("== Cross-car correlation (scripts/42; resampling the 1,198 ==")
w("   dual-coverage cell-hours)")
w(sprintf("  corr of per-car mean log(NOx): %.3f [%.3f, %.3f]",
          pt$cc_corr, ci(bs$cc_corr)[1], ci(bs$cc_corr)[2]))
w("")
w("== VERDICTS ==")
w(sprintf("  O3 composition gap excludes zero:  %s", excl0(bs$infl_o3)))
w(sprintf("  UFP composition gap excludes zero: %s", excl0(bs$infl_pn)))
w(sprintf("  Three pollutants mutually distinct (NOx-UFP and UFP-O3 both"))
w(sprintf("  exclude zero): %s",
          ifelse(excl0(bs$diff_nox_pn) == "YES" &
                 excl0(bs$diff_pn_o3) == "YES", "YES", "NO")))
close(out)
cat(readLines("output/tables/amt_bootstrap.txt"), sep = "\n")
