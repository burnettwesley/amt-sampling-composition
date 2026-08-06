# 00_generate_demo_data.R
# ---------------------------------------------------------------------------
# SYNTHETIC DEMONSTRATION DATA GENERATOR
#
# Builds a small stylized street network and ~120 collection days of two-car
# mobile-monitoring routes, then emits ~200k 1 Hz-style readings whose columns
# mimic the Aclima release schema. The output exists ONLY so that the paper's
# sampling-composition mechanism (revisit skew -> blind-join exposure
# inflation) can be demonstrated end-to-end by anyone, with no data agreement.
#
# COMPLIANCE RULE (absolute): this script is parameterized ONLY by summary
# statistics published in the manuscript (Burnett & Sheldon, "Sampling-
# composition bias in urban mobile air-quality monitoring", AMT). It reads
# NO Aclima file, no file under data/raw/, and no derived product of the
# release. Every number below is either (a) a published statistic, cited in
# the comment beside it, or (b) an arbitrary round convenience value,
# documented as such.
#
# Published statistics used (Oakland unless noted):
#   - Median revisit days per 100 m cell by class: 65 / 33 / 26 / 16
#     (freeway / arterial / collector / local), out of the campaign's
#     collection days.
#   - Mean log(NOx) by class: 4.26 / 3.42 / 3.23 / 2.69
#     (freeway / arterial / collector / local).
#   - Mean O3 (ppb, San Francisco release) by class, inverted by titration:
#     22.6 / 25.4 / 28.2 / 31.3 (freeway / arterial / collector / local).
#   - Particle number: shallow positive cross-class gradient (the paper's
#     composition artifact is +0.11 log points against +1.05 for NOx, i.e.
#     roughly one-tenth the NOx gradient).
#   - Freeway readings present on ~100 percent of collection days.
#
# Convenience values (NOT from any data file; round, documented choices):
#   - within-class SD of log(NOx): 0.80        (round value)
#   - within-class SD of O3:       6 ppb       (round value)
#   - within-class SD of log(PN):  0.50        (round value)
#   - diurnal amplitude of log(NOx): 0.25, cosine peaking 08:00 (qualitative)
#   - grid size, spacing, speeds, points-per-pass, NO2/NO split, PN bin
#     shares: stylized, chosen for a compact offline demo.
#
# Outputs:
#   demo_data/demo_readings.csv.gz  (~200k rows; first line is a disclaimer)
#   demo_data/demo_network.csv      (synthetic segments with road class)
#   demo_data/README.md             (provenance + compliance statement)
#
# Run from the repository root:  Rscript scripts/00_generate_demo_data.R
# Runs offline, base R only. Fixed seed; fully reproducible.
# ---------------------------------------------------------------------------

set.seed(20260806)

t_start <- Sys.time()
dir.create("demo_data", showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. Stylized street network: ~30x30 grid + one freeway spine
## ---------------------------------------------------------------------------
## Synthetic coordinates in a fake UTM-like box (meters). The box is placed,
## for the lon/lat columns only, over open water in Monterey Bay so that the
## coordinates cannot be mistaken for any real street network.

nx <- 30L; ny <- 30L          # grid nodes per side
spacing <- 100                 # meters between grid lines (one 100 m cell per segment)
x_orig <- 550000; y_orig <- 4180000   # fake UTM-like origin

grid_x <- x_orig + (0:(nx - 1L)) * spacing   # vertical-line x positions
grid_y <- y_orig + (0:(ny - 1L)) * spacing   # horizontal-line y positions

## Class layout (stylized):
##  - one freeway spine: a separate vertical line threaded between grid
##    columns 15 and 16 (offset so it is its own corridor)
##  - arterials: horizontal rows 8 and 23, vertical column 8
##  - collectors: horizontal rows 4 and 27, vertical columns 15 and 22
##  - all remaining grid lines are local streets
spine_x  <- x_orig + 14.5 * spacing + 37    # freeway spine x (between cols 15/16)
art_rows <- c(8L, 23L); art_cols <- c(8L)
col_rows <- c(4L, 27L); col_cols <- c(15L, 22L)

line_class <- function(idx, art_idx, col_idx) {
  ifelse(idx %in% art_idx, "arterial", ifelse(idx %in% col_idx, "collector", "local"))
}

## Build segment table: each grid-line is chopped into 29 segments of 100 m.
segs <- list()
## vertical grid lines
for (i in seq_len(nx)) {
  cls <- line_class(i, art_cols, col_cols)
  segs[[length(segs) + 1L]] <- data.frame(
    orientation = "v", line_coord = grid_x[i], seg_idx = 1:(ny - 1L),
    x0 = grid_x[i], y0 = grid_y[-ny], x1 = grid_x[i], y1 = grid_y[-1L],
    class = cls)
}
## horizontal grid lines
for (j in seq_len(ny)) {
  cls <- line_class(j, art_rows, col_rows)
  segs[[length(segs) + 1L]] <- data.frame(
    orientation = "h", line_coord = grid_y[j], seg_idx = 1:(nx - 1L),
    x0 = grid_x[-nx], y0 = grid_y[j], x1 = grid_x[-1L], y1 = grid_y[j],
    class = cls)
}
## freeway spine (its own vertical line)
segs[[length(segs) + 1L]] <- data.frame(
  orientation = "v", line_coord = spine_x, seg_idx = 1:(ny - 1L),
  x0 = spine_x, y0 = grid_y[-ny], x1 = spine_x, y1 = grid_y[-1L],
  class = "freeway")

net <- do.call(rbind, segs)
net$seg_id <- seq_len(nrow(net))
net$class <- factor(net$class, levels = c("freeway", "arterial", "collector", "local"))

cat(sprintf("Network: %d segments (%s)\n", nrow(net),
            paste(sprintf("%s %d", levels(net$class), as.integer(table(net$class))),
                  collapse = ", ")))

## ---------------------------------------------------------------------------
## 2. 120 collection days of two-car routes (revisit-skew calibration)
## ---------------------------------------------------------------------------
## Marginal per-day visit probabilities are chosen so that the MEDIAN
## distinct-visit-day count per segment matches the published Oakland medians
## 65 / 33 / 26 / 16 out of the campaign's days (here scaled to 120 synthetic
## collection days): p = published_median / 120.
##
## Routing structure (qualitative, per the published campaign description):
## drivers received DAILY POLYGON ASSIGNMENTS, so minor-street coverage is
## concentrated in one neighborhood at a time while the connecting drives put
## the fleet on the major-road skeleton citywide every day. We mimic this
## with a 3x3 grid of neighborhood zones: local and collector segments are
## visited only when their zone is active that day (zone-active probability
## a, within-zone visit probability q, with a*q = published_median/120), while
## freeway and arterial segments draw citywide every day. The freeway spine
## is guaranteed present every day (published: freeway readings on ~100
## percent of collection days).

n_days <- 120L
## 120 consecutive weekdays starting 2016-01-04 (synthetic calendar)
all_days <- seq(as.Date("2016-01-04"), by = "day", length.out = 200)
coll_days <- all_days[!format(all_days, "%u") %in% c("6", "7")][1:n_days]

p_visit <- c(freeway = 65 / 120,    # published median 65 revisit days
             arterial = 33 / 120,   # published median 33
             collector = 26 / 120,  # published median 26
             local = 16 / 120)      # published median 16 (~15% of days)

## neighborhood zones (3x3 over the grid, by segment midpoint)
extent <- (nx - 1L) * spacing
mid_x <- (net$x0 + net$x1) / 2; mid_y <- (net$y0 + net$y1) / 2
zc <- pmin(floor((mid_x - x_orig) / (extent / 3)), 2)
zr <- pmin(floor((mid_y - y_orig) / (extent / 3)), 2)
net$zone <- 1L + zc + 3L * zr

## zone-active probabilities (stylized routing intensities); within-zone visit
## probabilities follow from the calibration a * q = published_median / 120
a_zone <- c(collector = 0.30, local = 0.20)
q_zone <- c(collector = unname(p_visit["collector"] / a_zone["collector"]),
            local     = unname(p_visit["local"]     / a_zone["local"]))

visit <- matrix(FALSE, nrow(net), n_days)
is_major <- net$class %in% c("freeway", "arterial")
pv_major <- p_visit[as.character(net$class[is_major])]
visit[is_major, ] <- matrix(runif(sum(is_major) * n_days) < pv_major, sum(is_major))
for (cl in c("collector", "local")) {
  A <- matrix(runif(9L * n_days) < a_zone[cl], 9L, n_days)   # zone x day activity
  rows <- which(net$class == cl)
  active <- A[net$zone[rows], , drop = FALSE]                # active-zone indicator
  visit[rows, ] <- active & (matrix(runif(length(rows) * n_days),
                                    length(rows)) < q_zone[cl])
}

## guarantee freeway presence every day (probability of a miss is ~1e-10,
## but make it structural, as in the published campaign)
fw_ids <- which(net$class == "freeway")
for (d in which(colSums(visit[fw_ids, , drop = FALSE]) == 0L))
  visit[sample(fw_ids, 1L), d] <- TRUE

## achieved revisit medians (distinct visit days per 100 m segment, by class)
revisit <- rowSums(visit)
ach_med <- tapply(revisit, net$class, median)
cat("Achieved revisit-day medians (targets 65/33/26/16):\n")
print(round(ach_med, 1))
fw_day_share <- mean(colSums(visit[fw_ids, , drop = FALSE]) > 0L)
cat(sprintf("Share of days with freeway coverage: %.1f%% (published: 100%%)\n",
            100 * fw_day_share))

## ---------------------------------------------------------------------------
## 3. Emit 1 Hz-style readings along each pass
## ---------------------------------------------------------------------------
## Points per pass = segment length / speed at 1 Hz (stylized speeds; the
## local-street value is chosen high so the total lands near ~200k rows
## while the 15%-of-days local visit rate is preserved).
speed_ms <- c(freeway = 26, arterial = 16, collector = 14, local = 16)  # m/s, stylized
n_pts_cl <- pmax(round(spacing / speed_ms), 2L)                         # readings per pass

pass_idx <- which(visit, arr.ind = TRUE)                # (segment, day) pairs
passes <- data.frame(seg = pass_idx[, 1L], day = pass_idx[, 2L])
passes <- passes[order(passes$day, passes$seg), ]
passes$class <- net$class[passes$seg]
passes$n_pts <- n_pts_cl[as.character(passes$class)]
passes$car <- sample(c("Car_A", "Car_B"), nrow(passes), replace = TRUE)
## pass start time: uniform within a 09:00-17:00 working window (weekday
## working-hours collection, as in the published campaign description)
passes$start_s <- 9L * 3600L + floor(runif(nrow(passes)) * (8L * 3600L - 60L))

n_row <- sum(passes$n_pts)
cat(sprintf("Passes: %d; readings to emit: %d\n", nrow(passes), n_row))

## expand passes to points
rep_id <- rep(seq_len(nrow(passes)), passes$n_pts)
k      <- sequence(passes$n_pts)                        # 1..n_pts within pass
frac   <- (k - 0.5) / passes$n_pts[rep_id]              # position along segment
seg    <- passes$seg[rep_id]
x <- net$x0[seg] + frac * (net$x1[seg] - net$x0[seg]) + rnorm(n_row, 0, 2)  # 2 m GPS jitter
y <- net$y0[seg] + frac * (net$y1[seg] - net$y0[seg]) + rnorm(n_row, 0, 2)
day_i  <- passes$day[rep_id]
tsec   <- passes$start_s[rep_id] + (k - 1L)             # 1 Hz within pass
cls    <- as.character(net$class[seg])
hour   <- tsec %/% 3600

## --- pollutant draws (published class means; round convenience SDs) --------
mu_lognox <- c(freeway = 4.26, arterial = 3.42, collector = 3.23, local = 2.69)  # published
mu_o3     <- c(freeway = 22.6, arterial = 25.4, collector = 28.2, local = 31.3)  # published
## PN: shallow positive gradient, roughly one-tenth of the 1.57 log-point
## NOx freeway-local gap (published artifact ratio +0.11 vs +1.05), so a
## freeway-local log(PN) gap of 0.16 around an arbitrary base of 9.20:
mu_logpn  <- c(freeway = 9.36, arterial = 9.31, collector = 9.25, local = 9.20)

## diurnal term: qualitative morning-peak cosine, amplitude 0.25 (convenience),
## de-meaned over the 09-17 sampling window so class means stay on target
diurnal_raw <- function(h) 0.25 * cos(2 * pi * (h - 8) / 24)
diur <- diurnal_raw(hour) - mean(diurnal_raw(9:16))

log_nox <- mu_lognox[cls] + diur + rnorm(n_row, 0, 0.80)   # SD 0.80: round value
nox <- exp(log_nox)
r_no2 <- plogis(rnorm(n_row, qlogis(0.55), 0.30))          # NO2 share ~0.55, qualitative
NO2 <- nox * r_no2
NO  <- nox * (1 - r_no2)

O3 <- pmax(rnorm(n_row, mu_o3[cls] - 4 * diur, 6), 0.5)    # SD 6 ppb: round value;
                                                           # weak anti-phase diurnal (qualitative)
pn_tot <- exp(mu_logpn[cls] + rnorm(n_row, 0, 0.50))       # SD 0.50: round value
shares <- c(0.42, 0.26, 0.16, 0.10, 0.06)                  # stylized bin shares
sh_noise <- matrix(exp(rnorm(n_row * 5L, 0, 0.20)), ncol = 5L)
sh <- sweep(sh_noise, 2, shares, `*`); sh <- sh / rowSums(sh)
PN <- round(sh * pn_tot)

car_speed <- pmax(round(speed_ms[cls] + rnorm(n_row, 0, 2), 1), 1)

## --- fake lon/lat (equirectangular around a point in open water) -----------
lat0 <- 36.90; lon0 <- -121.90                              # Monterey Bay (water)
m_per_deg_lat <- 111132
m_per_deg_lon <- 111320 * cos(lat0 * pi / 180)
Latitude  <- round(lat0 + (y - y_orig) / m_per_deg_lat, 6)
Longitude <- round(lon0 + (x - x_orig) / m_per_deg_lon, 6)

Date_Time <- paste(format(coll_days[day_i], "%Y-%m-%d"),
                   sprintf("%02d:%02d:%02d", tsec %/% 3600, (tsec %% 3600) %/% 60, tsec %% 60))

readings <- data.frame(
  Date_Time = Date_Time,
  Latitude = Latitude, Longitude = Longitude,
  NO2 = round(NO2, 3), NO = round(NO, 3), O3 = round(O3, 2),
  PN1 = PN[, 1], PN2 = PN[, 2], PN3 = PN[, 3], PN4 = PN[, 4], PN5 = PN[, 5],
  Car_Identifier = passes$car[rep_id], Car_Speed = car_speed)
readings <- readings[order(readings$Date_Time, readings$Car_Identifier), ]

## achieved class means (sanity report; not read from anywhere)
ach_mu <- tapply(log_nox, cls, mean)
cat("Achieved mean log(NOx) by class (targets 4.26/3.42/3.23/2.69):\n")
print(round(ach_mu[c("freeway", "arterial", "collector", "local")], 2))

## ---------------------------------------------------------------------------
## 4. Write outputs
## ---------------------------------------------------------------------------
disclaimer <- "# SYNTHETIC DEMONSTRATION DATA - contains no Aclima records; structural properties only; not for scientific inference"

con <- gzfile("demo_data/demo_readings.csv.gz", "w")
writeLines(disclaimer, con)
write.csv(readings, con, row.names = FALSE, quote = FALSE)
close(con)

write.csv(net[, c("seg_id", "orientation", "line_coord", "seg_idx",
                  "x0", "y0", "x1", "y1", "class")],
          "demo_data/demo_network.csv", row.names = FALSE, quote = FALSE)

readme <- sprintf(
"# Synthetic demonstration data

**SYNTHETIC DEMONSTRATION DATA — contains no Aclima records; structural
properties only; not for scientific inference.**

## Provenance

Every row in `demo_readings.csv.gz` was generated by
`scripts/00_generate_demo_data.R` (fixed seed 20260806) on a stylized 30x30
street grid with one synthetic freeway spine, placed over open water in
Monterey Bay precisely so the coordinates cannot be mistaken for a real
network. No Aclima file, no file under `data/raw/`, and no derived product of
the Aclima/EDF release was read at any point in the generation.

## Published-stats-only rule (compliance basis)

The generator is parameterized ONLY by summary statistics published in the
manuscript (Burnett & Sheldon, *Sampling-composition bias in urban mobile
air-quality monitoring*, AMT):

- median revisit days per 100 m cell by road class (65/33/26/16,
  freeway/arterial/collector/local, Oakland), used as marginal per-day visit
  probabilities over 120 synthetic collection days (minor streets are visited
  through daily neighborhood-zone assignments, mimicking the published
  description of the campaign's daily polygon routing);
- mean log(NOx) by class (4.26/3.42/3.23/2.69, Oakland);
- mean O3 by class (22.6/25.4/28.2/31.3 ppb, San Francisco release; inverted
  ordering per NO titration);
- a shallow positive particle-number gradient (~one-tenth of the NOx
  gradient, per the published +0.11 vs +1.05 artifact ratio);
- freeway coverage on ~100%% of collection days.

All other values (within-class SDs 0.80 / 6 ppb / 0.50, diurnal amplitude
0.25, grid geometry, speeds, NO2/NO split, PN bin shares) are round,
documented convenience choices carrying no information about any individual
Aclima record. Aggregate statistics published in a manuscript are not subject
to the data agreement's redistribution restriction; no record-level data are
reproduced or reproducible from these files.

## Achieved calibration (this build)

- revisit-day medians by class: %s (targets 65/33/26/16)
- mean log(NOx) by class: %s (targets 4.26/3.42/3.23/2.69)
- readings: %d rows over %d collection days, 2 cars

## Files

- `demo_readings.csv.gz` — ~200k 1 Hz-style readings mimicking the release
  schema (`Date_Time, Latitude, Longitude, NO2, NO, O3, PN1..PN5,
  Car_Identifier, Car_Speed`); first line is the disclaimer above.
- `demo_network.csv` — the synthetic segment network with road classes.

Run `Rscript scripts/99_demo_run.R` for the end-to-end demonstration.
",
  paste(round(ach_med[c("freeway", "arterial", "collector", "local")], 0), collapse = "/"),
  paste(round(ach_mu[c("freeway", "arterial", "collector", "local")], 2), collapse = "/"),
  nrow(readings), n_days)
writeLines(readme, "demo_data/README.md")

sz <- file.size("demo_data/demo_readings.csv.gz") / 1e6
cat(sprintf("\nWrote demo_data/demo_readings.csv.gz (%.1f MB), demo_network.csv, README.md\n", sz))
cat(sprintf("Elapsed: %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
