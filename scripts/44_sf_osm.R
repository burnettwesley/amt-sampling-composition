# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 44_sf_osm.R
# ===========
# SF replication, step 1: road network and reading classification.
# Downloads the San Francisco OSM road network (Overpass; same class mapping
# as scripts/33), classifies every SF reading in the California Aclima
# release (Jun 2016 - Sep 2017; file in the companion repo) by nearest
# segment, and caches both.
#
# OUTPUT: data/processed/sf_osm_roads.gpkg
#         data/processed/sf_reading_road_class.csv.gz  (row-indexed to the
#           SF-bbox subset defined here; the subset definition is re-derivable
#           from the filters below)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(lubridate); library(jsonlite)
})

UTM10N <- 32610
CA_FILE <- Sys.getenv("AQ_CA_FILE", "data/raw/aclima/California_201605_201709_GoogleAclimaAQ.txt")
SF <- c(lat_min = 37.708, lat_max = 37.833, lon_min = -122.515, lon_max = -122.355)

osm_json    <- "data/raw/osm/sf_roads_osm.json"
roads_cache <- "data/processed/sf_osm_roads.gpkg"
class_out   <- "data/processed/sf_reading_road_class.csv.gz"

# -- 1. Network --------------------------------------------------------------
if (!file.exists(osm_json)) {
  cat("[1/3] Downloading SF OSM network via Overpass ...\n")
  q <- paste0(
    '[out:json][timeout:240];way["highway"~"^(motorway|motorway_link|trunk|',
    'trunk_link|primary|primary_link|secondary|secondary_link|tertiary|',
    'tertiary_link|residential|unclassified|living_street|service)$"]',
    '(37.69,-122.53,37.84,-122.34);out geom;')
  ok <- system2("curl",
    c("-s", "--max-time", "300",
      "-A", shQuote("oakland-airquality-research/1.0 (burnettwesley@gmail.com)"),
      "-o", shQuote(osm_json), "--data-urlencode", shQuote(paste0("data=", q)),
      "https://overpass-api.de/api/interpreter"))
  stopifnot(ok == 0, file.size(osm_json) > 1e6)
} else cat("[1/3] Using cached SF Overpass JSON.\n")

if (file.exists(roads_cache)) {
  roads <- st_read(roads_cache, quiet = TRUE)
  cat("    cached network loaded.\n")
} else {
  osm <- fromJSON(osm_json, simplifyVector = FALSE)
  hw <- vapply(osm$elements, function(e) e$tags$highway %||% NA_character_,
               character(1))
  geoms <- lapply(osm$elements, function(e) {
    m <- t(vapply(e$geometry, function(g) c(g$lon, g$lat), numeric(2)))
    st_linestring(m)
  })
  roads <- st_sf(highway = hw, geometry = st_sfc(geoms, crs = 4326)) %>%
    filter(!is.na(highway)) %>%
    mutate(road_class = case_when(
      highway %in% c("motorway", "motorway_link")           ~ "freeway",
      highway %in% c("trunk", "trunk_link",
                     "primary", "primary_link")             ~ "arterial",
      highway %in% c("secondary", "secondary_link",
                     "tertiary", "tertiary_link")           ~ "collector",
      TRUE                                                  ~ "local")) %>%
    st_transform(UTM10N)
  st_write(roads, roads_cache, driver = "GPKG", quiet = TRUE, delete_dsn = TRUE)
}
cat(sprintf("    %s segments (%s freeway, %s arterial, %s collector, %s local)\n",
            format(nrow(roads), big.mark = ","),
            sum(roads$road_class == "freeway"),
            sum(roads$road_class == "arterial"),
            sum(roads$road_class == "collector"),
            sum(roads$road_class == "local")))

# network road-km by class (for readings-vs-network comparison)
roads$len_km <- as.numeric(st_length(roads)) / 1000
km <- roads %>% st_drop_geometry() %>% group_by(road_class) %>%
  summarise(km = sum(len_km), .groups = "drop") %>%
  mutate(pct = 100 * km / sum(km))
cat("    network road-km shares:\n")
for (i in seq_len(nrow(km)))
  cat(sprintf("      %-9s %6.0f km (%4.1f%%)\n", km$road_class[i],
              km$km[i], km$pct[i]))
write_csv(km, "data/processed/sf_network_km.csv")

# -- 2. SF readings ----------------------------------------------------------
cat("[2/3] Loading SF readings ...\n")
sf_pts <- read_csv(CA_FILE,
                   col_select = c(Date_Time, Latitude, Longitude, NO2, NO,
                                  Car_Identifier),
                   show_col_types = FALSE) %>%
  mutate(ts  = ymd_hms(Date_Time, tz = "UTC"),
         NO2 = suppressWarnings(as.numeric(NO2)),
         NO  = suppressWarnings(as.numeric(NO))) %>%
  filter(!is.na(Latitude), !is.na(Longitude), !is.na(ts),
         Latitude >= SF["lat_min"], Latitude <= SF["lat_max"],
         Longitude >= SF["lon_min"], Longitude <= SF["lon_max"]) %>%
  mutate(reading_uid = row_number())
cat(sprintf("    %s SF readings.\n", format(nrow(sf_pts), big.mark = ",")))

# -- 3. Classify -------------------------------------------------------------
cat("[3/3] Nearest-segment classification (chunked) ...\n")
CHUNK <- 250000L
n <- nrow(sf_pts); cls <- character(n); dst <- numeric(n)
pts_utm <- sf_pts %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(UTM10N)
for (k in seq_len(ceiling(n / CHUNK))) {
  i1 <- (k - 1L) * CHUNK + 1L; i2 <- min(k * CHUNK, n)
  cat(sprintf("    chunk %d / %d ...\n", k, ceiling(n / CHUNK)))
  sub <- pts_utm[i1:i2, ]
  idx <- st_nearest_feature(sub, roads)
  cls[i1:i2] <- roads$road_class[idx]
  dst[i1:i2] <- as.numeric(st_distance(sub, roads[idx, ], by_element = TRUE))
}
out <- tibble(reading_uid = sf_pts$reading_uid, road_class = cls,
              dist_road_m = round(dst, 1))
write_csv(out, class_out)
cat(sprintf("Saved %s (median snap %.1f m; freeway share %.1f%%)\n",
            class_out, median(out$dist_road_m),
            100 * mean(out$road_class == "freeway")))
