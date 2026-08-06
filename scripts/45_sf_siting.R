# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 45_sf_siting.R
# ==============
# SF replication, step 2: siting-composition diagnostics and the
# synthetic-outcome matching experiment, run identically for San Francisco
# (CA release, Jun 2016 - Sep 2017, route-free protocol) and Oakland
# (dedicated release, Jun 2015 - May 2016, polygon protocol).
#
# DIAGNOSTICS (readings only): road-class shares of readings vs network
# road-km; revisit frequency (~100m cells, distinct visit days) by class;
# mean log(NOx) by class; share of driving days with any freeway reading.
#
# SYNTHETIC-OUTCOME EXPERIMENT (the crash-free artifact demonstration):
# sample points uniformly at random along NON-freeway street segments;
# assign each a random timestamp (uniform driving day; hour drawn from the
# Oakland crash hour-of-day distribution, identical for both cities); match
# each point-time to its nearest reading under the paper's combined metric
# (2 km / +/-72 h, d = sqrt((m/500)^2 + (hr/6)^2)), blind and
# class-constrained; report the blind freeway-assignment share and the
# blind-minus-constrained gap in mean log(NOx). Oakland's blind
# freeway-assignment share validates the method against the real
# crash-based figure (59.1% of non-highway crashes received freeway
# readings). Seed fixed for reproducibility.
#
# OUTPUT: output/tables/sf_siting.txt
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(lubridate)
})
set.seed(20260720)

UTM10N <- 32610
N_SYN  <- 2000
CA_FILE <- Sys.getenv("AQ_CA_FILE", "data/raw/aclima/California_201605_201709_GoogleAclimaAQ.txt")
SFBB <- c(lat_min = 37.708, lat_max = 37.833, lon_min = -122.515, lon_max = -122.355)

# Oakland crash hour-of-day distribution (UTC), for outcome-time draws
cr_hours <- read_csv("data/processed/nearest_matched.csv",
                     show_col_types = FALSE) %>%
  mutate(h = hour(as.POSIXct(timestamp, tz = "UTC"))) %>% pull(h)

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
    km <- read_csv("data/processed/sf_network_km.csv", show_col_types = FALSE)
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
    roads$len_km <- as.numeric(st_length(roads)) / 1000
    km <- roads %>% st_drop_geometry() %>% group_by(road_class) %>%
      summarise(km = sum(len_km), .groups = "drop") %>%
      mutate(pct = 100 * km / sum(km))
  }
  xy <- pts %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(UTM10N) %>% st_coordinates()
  pts$x <- xy[, 1]; pts$y <- xy[, 2]
  pts <- pts %>%
    mutate(nox = NO2 + NO, date = as.Date(ts),
           cell = paste(floor(Longitude / 0.0012), floor(Latitude / 0.0009)))
  list(pts = pts, roads = roads, km = km)
}

diagnostics <- function(city, d, out) {
  pts <- d$pts
  cat(sprintf("== %s: composition diagnostics ==\n", city), file = out)
  sh <- pts %>% count(road_class) %>% mutate(pct = 100 * n / sum(n))
  kmt <- d$km %>% select(road_class, km_pct = pct)
  tab <- sh %>% left_join(kmt, by = "road_class")
  for (i in seq_len(nrow(tab)))
    cat(sprintf("  %-9s readings %5.1f%%  | network km %5.1f%%\n",
                tab$road_class[i], tab$pct[i], tab$km_pct[i]), file = out)
  rv <- pts %>% group_by(cell) %>%
    summarise(n_days = n_distinct(date),
              class = names(sort(table(road_class), decreasing = TRUE))[1],
              .groups = "drop") %>%
    group_by(class) %>% summarise(med = median(n_days), .groups = "drop")
  cat("  revisit days per ~100m cell (median): ", file = out)
  cat(paste(sprintf("%s %d", rv$class, rv$med), collapse = ", "), "\n",
      file = out)
  fw_days <- 100 * n_distinct(pts$date[pts$road_class == "freeway"]) /
    n_distinct(pts$date)
  cat(sprintf("  driving days with any freeway reading: %.1f%%\n", fw_days),
      file = out)
  mn <- pts %>% filter(!is.na(nox), nox > 0) %>%
    group_by(road_class) %>%
    summarise(m = mean(log(nox)), .groups = "drop")
  cat("  mean log(NOx) by class: ", file = out)
  cat(paste(sprintf("%s %.2f", mn$road_class, mn$m), collapse = ", "), "\n\n",
      file = out)
}

synthetic <- function(city, d, out) {
  pts <- d$pts %>% filter(!is.na(nox), nox > 0)
  roads <- d$roads
  # sample synthetic outcome sites along non-freeway segments, length-weighted
  nonf <- roads %>% filter(road_class != "freeway")
  samp <- st_sample(nonf, size = ceiling(N_SYN * 1.3), type = "random")
  samp <- st_cast(samp, "POINT")   # st_sample on lines can return MULTIPOINT
  sxy_full <- st_coordinates(samp)[, 1:2, drop = FALSE]
  good <- which(is.finite(sxy_full[, 1]) & is.finite(sxy_full[, 2]))
  stopifnot(length(good) >= N_SYN)
  good <- good[seq_len(N_SYN)]
  samp <- samp[good]
  sxy <- sxy_full[good, , drop = FALSE]
  # class of each synthetic site = nearest segment's class
  sidx <- st_nearest_feature(st_sf(geometry = samp), roads)
  s_class <- roads$road_class[sidx]
  stopifnot(length(s_class) == N_SYN, nrow(sxy) == N_SYN)
  # outcome times: uniform driving day, hour from Oakland crash distribution
  days <- sort(unique(pts$date))
  s_date <- sample(days, N_SYN, replace = TRUE)
  s_hour <- sample(cr_hours, N_SYN, replace = TRUE)
  s_ts <- as.POSIXct(paste(s_date), tz = "UTC") + s_hour * 3600 +
    runif(N_SYN, 0, 3600)

  match_one <- function(i, constrained) {
    t0 <- s_ts[i]
    cand <- pts %>%
      filter(ts >= t0 - 72 * 3600, ts <= t0 + 72 * 3600)
    if (constrained) cand <- cand %>% filter(road_class == s_class[i])
    if (nrow(cand) == 0) return(c(NA, NA))
    dx <- cand$x - sxy[i, 1]; dy <- cand$y - sxy[i, 2]
    sp <- sqrt(dx^2 + dy^2)
    keep <- sp <= 2000
    if (!any(keep)) return(c(NA, NA))
    cand <- cand[keep, ]; sp <- sp[keep]
    tmp <- abs(as.numeric(difftime(cand$ts, t0, units = "hours")))
    dd <- sqrt((sp / 500)^2 + (tmp / 6)^2)
    j <- which.min(dd)
    c(match(cand$road_class[j],
            c("local", "collector", "arterial", "freeway")),
      log(cand$nox[j]))
  }
  res_b <- t(vapply(seq_len(N_SYN), match_one, numeric(2), constrained = FALSE))
  res_c <- t(vapply(seq_len(N_SYN), match_one, numeric(2), constrained = TRUE))
  cls_lab <- c("local", "collector", "arterial", "freeway")
  ok_b <- !is.na(res_b[, 1]); ok_c <- !is.na(res_c[, 1])
  cat(sprintf("== %s: synthetic-outcome matching experiment (N = %d non-freeway sites) ==\n",
              city, N_SYN), file = out)
  cat(sprintf("  match rate: blind %.1f%%, class-constrained %.1f%%\n",
              100 * mean(ok_b), 100 * mean(ok_c)), file = out)
  cat(sprintf("  blind: share assigned a FREEWAY reading: %.1f%%\n",
              100 * mean(res_b[ok_b, 1] == 4)), file = out)
  own <- cls_lab[res_b[, 1]] == s_class
  cat(sprintf("  blind: share assigned a reading of the site's own class: %.1f%%\n",
              100 * mean(own[ok_b])), file = out)
  # Selection-free composition inflation: blind-assigned class mix valued at
  # citywide class means, minus the sites' own-class mix at the same means.
  cm <- pts %>% group_by(road_class) %>%
    summarise(m = mean(log(nox)), .groups = "drop")
  mvec <- setNames(cm$m, cm$road_class)
  infl <- mean(mvec[cls_lab[res_b[ok_b, 1]]]) - mean(mvec[s_class[ok_b]])
  cat(sprintf("  composition-implied inflation (blind class mix vs own-class mix, at citywide class means): %.2f log points\n",
              infl), file = out)
  both <- ok_b & ok_c
  gap <- mean(res_b[both, 2]) - mean(res_c[both, 2])
  cat(sprintf("  paired blind-vs-constrained gap (both-matched subset, n = %d): blind %.2f, constrained %.2f (gap %.2f log points)\n\n",
              sum(both), mean(res_b[both, 2]), mean(res_c[both, 2]), gap),
      file = out)
  invisible(NULL)
}

out <- file("output/tables/sf_siting.txt", open = "wt")
cat("SF siting-composition replication -- scripts/45 --",
    format(Sys.Date()), "\n", file = out)
cat(sprintf("Synthetic sites: %d per city, non-freeway segments, length-weighted;\n", N_SYN),
    file = out)
cat("times: uniform driving day x Oakland crash hour distribution; matching\n",
    file = out)
cat("rule identical to scripts/07 (2 km / +/-72 h combined metric).\n",
    file = out)
cat("Oakland validation target: 59.1% of real non-highway crashes received\n",
    file = out)
cat("freeway readings under blind matching.\n\n", file = out)

for (city in c("SF", "Oakland")) {
  d <- load_city(city)
  diagnostics(city, d, out)
  synthetic(city, d, out)
  rm(d); invisible(gc())
}
close(out)
cat(readLines("output/tables/sf_siting.txt"), sep = "\n")
