# 33_reading_road_class.R
# =======================
# Prep for road-class-aware matching (branch: roadclass-rematch).
# Classifies (a) every Aclima reading and (b) every study-period crash by the
# functional class of the nearest OSM road segment (freeway / arterial /
# collector / local; same mapping as scripts/32).
#
# Motivation (see SESSION_NOTES 2026-07-15): nearest-reading matching sited
# 59% of non-highway crashes on freeway readings because the fleet's reading
# supply is densest on freeways, mechanically inflating case exposure. The
# fix is to constrain matching to same-road-class readings; this script
# builds the class lookups that scripts/07 and /09 consume.
#
# Row-number contract: reading_uid = row_number() after script 07's cleaning
# filter (!is.na(Latitude/Longitude/timestamp)) on the raw Aclima file. Both
# 07 and 09 apply this exact filter before any other operation, so joining
# on reading_uid is exact in both. Crash numbering replicates script 07
# steps 3-4 verbatim and is VERIFIED against nearest_matched.csv.
#
# OUTPUT
# ------
#   data/processed/osm_roads.gpkg              (parsed/classified network)
#   data/processed/reading_road_class.csv.gz   (reading_uid, road_class, dist)
#   data/processed/crash_road_class.csv        (crash_id, road_class, dist,
#                                                STATE_HWY_IND cross-tab printed)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(lubridate)
  library(jsonlite)
})

UTM10N <- 32610

pollution_file <- "data/raw/aclima/Oakland_201505_201605_GoogleAclimaAQ.txt"
crashes_file   <- "data/raw/accidents/Crashes_combined.csv"
tract_cache    <- "data/processed/tracts.gpkg"
osm_file       <- "data/raw/osm/oakland_roads_osm.json"
roads_cache    <- "data/processed/osm_roads.gpkg"
reading_out    <- "data/processed/reading_road_class.csv.gz"
crash_out      <- "data/processed/crash_road_class.csv"

CHUNK <- 250000L


# == Step 1: road network =====================================================

if (file.exists(roads_cache)) {
  cat("[1/4] Loading cached road network ...\n")
  roads <- st_read(roads_cache, quiet = TRUE)
} else {
  cat("[1/4] Parsing OSM network ...\n")
  osm <- fromJSON(osm_file, simplifyVector = FALSE)
  hw <- vapply(osm$elements, function(e) e$tags$highway %||% NA_character_,
               character(1))
  geoms <- lapply(osm$elements, function(e) {
    m <- t(vapply(e$geometry, function(g) c(g$lon, g$lat), numeric(2)))
    st_linestring(m)
  })
  roads <- st_sf(highway = hw, geometry = st_sfc(geoms, crs = 4326)) %>%
    filter(!is.na(highway)) %>%
    mutate(road_class = case_when(
      highway %in% c("motorway", "motorway_link")            ~ "freeway",
      highway %in% c("trunk", "trunk_link",
                     "primary", "primary_link")              ~ "arterial",
      highway %in% c("secondary", "secondary_link",
                     "tertiary", "tertiary_link")            ~ "collector",
      TRUE                                                   ~ "local"
    )) %>%
    st_transform(UTM10N)
  st_write(roads, roads_cache, driver = "GPKG", quiet = TRUE, delete_dsn = TRUE)
}
cat(sprintf("    %s segments.\n", format(nrow(roads), big.mark = ",")))

classify_points <- function(x, y) {
  n <- length(x)
  cls  <- character(n)
  dst  <- numeric(n)
  n_chunks <- ceiling(n / CHUNK)
  for (k in seq_len(n_chunks)) {
    i1 <- (k - 1L) * CHUNK + 1L
    i2 <- min(k * CHUNK, n)
    if (n_chunks > 1) cat(sprintf("    chunk %d / %d ...\n", k, n_chunks))
    pts <- st_as_sf(data.frame(x = x[i1:i2], y = y[i1:i2]),
                    coords = c("x", "y"), crs = UTM10N)
    idx <- st_nearest_feature(pts, roads)
    cls[i1:i2] <- roads$road_class[idx]
    dst[i1:i2] <- as.numeric(st_distance(pts, roads[idx, ], by_element = TRUE))
  }
  list(class = cls, dist = dst)
}


# == Step 2: classify all readings ============================================

cat("[2/4] Classifying all Aclima readings (script 07 filter, reading_uid) ...\n")

pollution <- read_csv(pollution_file, show_col_types = FALSE) %>%
  mutate(timestamp = ymd_hms(Date_Time, tz = "UTC")) %>%
  filter(!is.na(Latitude), !is.na(Longitude), !is.na(timestamp))

cat(sprintf("    %s readings after cleaning filter.\n",
            format(nrow(pollution), big.mark = ",")))

poll_sf <- pollution %>%
  select(Latitude, Longitude) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(UTM10N)
pc <- st_coordinates(poll_sf)

res <- classify_points(pc[, 1], pc[, 2])

reading_class <- tibble(
  reading_uid = seq_len(nrow(pollution)),
  road_class  = res$class,
  dist_road_m = round(res$dist, 1)
)
write_csv(reading_class, reading_out)
cat(sprintf("    Saved %s (median snap %.1f m; freeway share %.1f%%).\n",
            reading_out, median(reading_class$dist_road_m),
            100 * mean(reading_class$road_class == "freeway")))


# == Step 3: classify crashes (script 07 steps 3-4, verbatim numbering) ======

cat("[3/4] Classifying crashes ...\n")

crashes <- read_csv(crashes_file, show_col_types = FALSE,
                    col_types = cols(.default = "c")) %>%
  mutate(
    COLLISION_TIME_int = suppressWarnings(as.integer(COLLISION_TIME)),
    hh        = COLLISION_TIME_int %/% 100L,
    mm        = COLLISION_TIME_int %%  100L,
    timestamp = ymd_hm(
      paste0(COLLISION_DATE, " ", sprintf("%02d", hh), ":", sprintf("%02d", mm)),
      tz = "UTC"
    ),
    lat = suppressWarnings(as.numeric(coalesce(LATITUDE,  POINT_Y))),
    lon = suppressWarnings(as.numeric(coalesce(LONGITUDE, POINT_X)))
  ) %>%
  filter(!is.na(lat), !is.na(lon), !is.na(timestamp)) %>%
  filter(
    timestamp >= ymd("2015-06-01", tz = "UTC"),
    timestamp <= ymd_hms("2016-05-31 23:59:59", tz = "UTC")
  ) %>%
  mutate(crash_id = row_number())

crash_sf <- crashes %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(UTM10N)
cc <- st_coordinates(crash_sf)

resc <- classify_points(cc[, 1], cc[, 2])

crash_class <- tibble(
  crash_id        = crashes$crash_id,
  crash_road_class = resc$class,
  crash_dist_road_m = round(resc$dist, 1),
  STATE_HWY_IND   = crashes$STATE_HWY_IND
)

# -- Verify crash_id numbering against the existing matched file -------------
chk <- read_csv("data/processed/nearest_matched.csv", show_col_types = FALSE) %>%
  select(crash_id, crash_x, crash_y)
chk_join <- chk %>%
  left_join(tibble(crash_id = crashes$crash_id,
                   x = cc[, 1], y = cc[, 2]), by = "crash_id")
max_dev <- max(abs(chk_join$crash_x - chk_join$x),
               abs(chk_join$crash_y - chk_join$y), na.rm = TRUE)
cat(sprintf("    crash_id verification vs nearest_matched.csv: max coord dev %.3f m (n=%s)\n",
            max_dev, format(nrow(chk_join), big.mark = ",")))
stopifnot(max_dev < 0.5)

write_csv(crash_class, crash_out)


# == Step 4: diagnostics ======================================================

cat("[4/4] Crash class distribution and STATE_HWY_IND cross-tab:\n")
print(crash_class %>% count(crash_road_class) %>% mutate(pct = round(100 * n / sum(n), 1)))
print(crash_class %>%
        count(STATE_HWY_IND, crash_road_class) %>%
        group_by(STATE_HWY_IND) %>% mutate(pct = round(100 * n / sum(n), 1)) %>%
        ungroup())
cat(sprintf("median crash snap distance: %.1f m (90th pctile %.1f m)\n",
            median(crash_class$crash_dist_road_m),
            quantile(crash_class$crash_dist_road_m, 0.9)))
cat("Done.\n")
