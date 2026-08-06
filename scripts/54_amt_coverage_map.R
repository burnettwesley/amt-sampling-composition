# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 54_amt_coverage_map.R
# =====================
# AMT Figure 1: two-panel map of sampling coverage colored by revisit
# frequency (distinct collection days per ~100 m cell), Oakland and San
# Francisco. Contains NO outcome data (bright line). Visualizes the
# revisit-skew claim of AMT section 3.1.
# OUTPUT: output/figures/fig_amt_coverage.png
# =============================================================================

suppressPackageStartupMessages({library(tidyverse); library(lubridate); library(sf)})
UTM10N <- 32610; CELL_M <- 100

cellize <- function(df) {
  pts <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326) %>% st_transform(UTM10N)
  xy <- st_coordinates(pts)
  df %>% mutate(cx = floor(xy[,1] / CELL_M), cy = floor(xy[,2] / CELL_M)) %>%
    group_by(cx, cy) %>% summarise(days = n_distinct(day), .groups = "drop") %>%
    mutate(x = cx * CELL_M + CELL_M/2, y = cy * CELL_M + CELL_M/2)
}
to_ll <- function(cells) {
  p <- st_as_sf(cells, coords = c("x", "y"), crs = UTM10N) %>% st_transform(4326)
  cbind(cells, st_coordinates(p) %>% as_tibble() %>% rename(lon = X, lat = Y))
}

oak <- read_csv("data/raw/aclima/Oakland_201505_201605_GoogleAclimaAQ.txt",
                col_select = c(Date_Time, Latitude, Longitude), show_col_types = FALSE) %>%
  filter(!is.na(Latitude)) %>%
  transmute(lat = Latitude, lon = Longitude,
            day = as.Date(ymd_hms(Date_Time, tz = "UTC"))) %>%
  cellize() %>% to_ll() %>% mutate(city = "Oakland (2015-2016)")

sfd <- read_csv(Sys.getenv("AQ_CA_FILE", "data/raw/aclima/California_201605_201709_GoogleAclimaAQ.txt"),
                col_select = c(Date_Time, Latitude, Longitude), show_col_types = FALSE) %>%
  filter(!is.na(Latitude), Latitude > 37.70, Latitude < 37.84,
         Longitude > -122.53, Longitude < -122.35) %>%
  transmute(lat = Latitude, lon = Longitude,
            day = as.Date(ymd_hms(Date_Time, tz = "UTC"))) %>%
  cellize() %>% to_ll() %>% mutate(city = "San Francisco (2016-2017)")

both <- bind_rows(oak, sfd) %>%
  mutate(days_c = cut(days, c(0, 5, 15, 30, 60, Inf),
                      labels = c("1-5", "6-15", "16-30", "31-60", "> 60")))

g <- ggplot(both, aes(lon, lat, color = days_c)) +
  geom_point(size = 0.5, shape = 15) +
  scale_color_viridis_d(name = "Distinct collection\ndays per 100 m cell",
                        option = "inferno", direction = -1, begin = 0.10, end = 0.85) +
  facet_wrap(~city, scales = "free") +
  coord_quickmap() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_minimal(base_size = 11) +
  theme(axis.title = element_blank(), panel.grid = element_blank(),
        legend.position = "right")

ggsave("output/figures/fig_amt_coverage.png", g, width = 11, height = 5.2, dpi = 300, bg = "white")
cat("cells:", nrow(oak), "Oakland,", nrow(sfd), "SF\n")
cat("share of cells with >60 days: Oakland",
    round(100*mean(oak$days > 60), 1), "% | SF", round(100*mean(sfd$days > 60), 1), "%\n")
