# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 68_distance_decay.R
# ===================
# Distance-decay of NOx away from the freeway network, for the AMT paper's
# cross-class gradient subsection (Sect. 3.2). Grounds the "steep gradient"
# premise in our own data: mean log(NOx) by distance-to-nearest-freeway bin
# (0-50, 50-100, 100-200, 200-400, 400-800, >800 m), NON-FREEWAY readings
# only, in both deployments.
#
# Machinery: the cached OSM networks from scripts/33 (Oakland) and /44 (SF)
# supply the classified segments; reading road classes come from the cached
# reading_road_class files (row-aligned to the script-07 / script-44 cleaning
# filters, verified by nrow stopifnot as in scripts/56). Distances to the
# freeway subnetwork are computed fresh here (the cached dist_road_m is
# distance to the nearest road of ANY class, not to the freeway).
#
# Reported per city: mean log(NOx) and n per bin; total decay (bin 1 minus
# last bin, log points); and the half-decay distance, i.e., the distance at
# which half the freeway increment (bin-1 mean minus >800 m far-field mean)
# has decayed, linearly interpolated across bin midpoints.
#
# OUTPUT: output/tables/distance_decay.txt
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(lubridate)
})

UTM10N  <- 32610
CA_FILE <- Sys.getenv("AQ_CA_FILE", "data/raw/aclima/California_201605_201709_GoogleAclimaAQ.txt")
SFBB    <- c(lat_min = 37.708, lat_max = 37.833,
             lon_min = -122.515, lon_max = -122.355)
CHUNK   <- 250000L

BREAKS <- c(0, 50, 100, 200, 400, 800, Inf)
LABELS <- c("0-50", "50-100", "100-200", "200-400", "400-800", ">800")
MIDS   <- c(25, 75, 150, 300, 600, NA)   # last bin open-ended; anchor at 1000
MIDS[6] <- 1000

# -- loading (filters verbatim from scripts/56) -------------------------------
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
  pts <- pts %>%
    mutate(nox = NO2 + NO) %>%
    filter(!is.na(nox), nox > 0, road_class != "freeway")
  list(pts = pts, freeway = roads %>% filter(road_class == "freeway"))
}

# -- distance of every reading to the freeway subnetwork (chunked) ------------
freeway_dist <- function(pts, freeway) {
  xy <- pts %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
    st_transform(UTM10N)
  n <- nrow(pts); dst <- numeric(n)
  for (k in seq_len(ceiling(n / CHUNK))) {
    i1 <- (k - 1L) * CHUNK + 1L; i2 <- min(k * CHUNK, n)
    cat(sprintf("    chunk %d / %d ...\n", k, ceiling(n / CHUNK)))
    sub <- xy[i1:i2, ]
    idx <- st_nearest_feature(sub, freeway)
    dst[i1:i2] <- as.numeric(st_distance(sub, freeway[idx, ],
                                         by_element = TRUE))
  }
  dst
}

decay_stats <- function(city) {
  cat(sprintf("== %s ==\n", city))
  d <- load_city(city)
  cat(sprintf("    %s complete-NOx non-freeway readings; %d freeway segments\n",
              format(nrow(d$pts), big.mark = ","), nrow(d$freeway)))
  d$pts$dist_fw <- freeway_dist(d$pts, d$freeway)
  tab <- d$pts %>%
    mutate(bin = cut(dist_fw, BREAKS, labels = LABELS, right = FALSE)) %>%
    group_by(bin) %>%
    summarise(mean_lnox = mean(log(nox)), n = n(), .groups = "drop") %>%
    arrange(bin)
  stopifnot(nrow(tab) == length(LABELS))
  total <- tab$mean_lnox[1] - tab$mean_lnox[nrow(tab)]
  # half-decay: distance where mean = far-field + 0.5 * increment,
  # interpolated linearly across bin midpoints
  target <- tab$mean_lnox[nrow(tab)] + 0.5 * total
  half <- approx(x = tab$mean_lnox, y = MIDS, xout = target, ties = mean)$y
  list(tab = tab, total = total, half = half)
}

res <- list(Oakland = decay_stats("Oakland"), SF = decay_stats("SF"))

# -- archive ------------------------------------------------------------------
out <- file("output/tables/distance_decay.txt", open = "wt")
w <- function(...) cat(..., "\n", sep = "", file = out)
w("NOx distance-decay from the freeway network -- scripts/68 -- ",
  format(Sys.Date()))
w("Mean log(NOx) by distance to the nearest OSM freeway segment (motorway/")
w("motorway_link; cached networks from scripts/33 and /44), NON-FREEWAY")
w("readings only, complete-NOx sample (NO2 + NO observed and positive;")
w("filters verbatim from scripts/56). Distances computed by snapping each")
w("reading to the freeway subnetwork in UTM 10N. Half-decay distance: where")
w("half the freeway increment (bin-1 mean minus >800 m far-field mean) has")
w("decayed, linearly interpolated across bin midpoints (25, 75, 150, 300,")
w("600, 1000 m anchor for the open bin).")
w("")
for (city in names(res)) {
  r <- res[[city]]
  w(sprintf("== %s ==", city))
  w(sprintf("  %-8s %10s %12s", "bin (m)", "mean lnox", "n"))
  for (i in seq_len(nrow(r$tab)))
    w(sprintf("  %-8s %10.3f %12s", as.character(r$tab$bin[i]),
              r$tab$mean_lnox[i], format(r$tab$n[i], big.mark = ",")))
  w(sprintf("  total decay (0-50 m minus >800 m): %.2f log points", r$total))
  w(sprintf("  half-decay distance: %.0f m", r$half))
  w("")
}
close(out)
cat(readLines("output/tables/distance_decay.txt"), sep = "\n")
