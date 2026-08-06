# 99_demo_run.R
# ---------------------------------------------------------------------------
# END-TO-END DEMONSTRATION on the synthetic dataset (no data agreement
# required). Mirrors the paper's pipeline at reduced scale:
#
#   1. classify every reading by the synthetic network's road classes;
#   2. run the pre-join diagnostics (readings share vs network share;
#      revisit-day medians by class);
#   3. synthetic-outcome experiment: 500 random non-freeway sites, matched
#      blind vs class-constrained under the paper's join metric
#      d = sqrt((ds/500 m)^2 + (dt/6 h)^2), caps 2 km and +/-72 h;
#   4. nonparametric bootstrap over sites (B = 199, percentile 95% CIs).
#
# Expected qualitative result (the paper's artifact): POSITIVE composition-
# implied log(NOx) inflation at non-freeway sites, an INVERTED (negative) O3
# gap, and a small positive particle-number gap.
#
# SYNTHETIC DEMONSTRATION DATA - contains no Aclima records; structural
# properties only; not for scientific inference.
#
# Run from the repository root, after scripts/00_generate_demo_data.R:
#   Rscript scripts/99_demo_run.R
# Offline, base R only; completes in well under 5 minutes.
# ---------------------------------------------------------------------------

t_start <- Sys.time()
set.seed(20260807)

## ---------------------------------------------------------------------------
## 0. Load synthetic data
## ---------------------------------------------------------------------------
net <- read.csv("demo_data/demo_network.csv", stringsAsFactors = FALSE)
con <- gzfile("demo_data/demo_readings.csv.gz")
rd  <- read.csv(con, skip = 1, stringsAsFactors = FALSE)   # skip disclaimer line
cat(sprintf("Loaded %d readings, %d network segments\n", nrow(rd), nrow(net)))

## recover planar coordinates (same fake equirectangular constants as 00_)
lat0 <- 36.90; lon0 <- -121.90
x_orig <- 550000; y_orig <- 4180000
m_per_deg_lat <- 111132
m_per_deg_lon <- 111320 * cos(lat0 * pi / 180)
rd$x <- x_orig + (rd$Longitude - lon0) * m_per_deg_lon
rd$y <- y_orig + (rd$Latitude  - lat0) * m_per_deg_lat
rd$t <- as.numeric(as.POSIXct(rd$Date_Time, tz = "UTC", format = "%Y-%m-%d %H:%M:%S"))
rd$date <- substr(rd$Date_Time, 1, 10)

## ---------------------------------------------------------------------------
## 1. Classify readings by the synthetic network's classes (snap to nearest line)
## ---------------------------------------------------------------------------
## The network is a set of axis-parallel lines chopped into 100 m segments,
## so the nearest segment is found exactly: nearest vertical line vs nearest
## horizontal line, then the segment index from the position along the line.
vx <- sort(unique(net$line_coord[net$orientation == "v"]))
hy <- sort(unique(net$line_coord[net$orientation == "h"]))
spacing <- diff(sort(unique(net$y0[net$orientation == "v"])))[1]
y_min <- min(net$y0[net$orientation == "v"])
x_min <- min(net$x0[net$orientation == "h"])
n_seg_line <- max(net$seg_idx)

nearest <- function(q, grid) {           # nearest grid value for each query
  i <- findInterval(q, grid, all.inside = TRUE)
  lo <- grid[i]; hi <- grid[pmin(i + 1L, length(grid))]
  ifelse(abs(q - lo) <= abs(q - hi), lo, hi)
}
nvx <- nearest(rd$x, vx); dvx <- abs(rd$x - nvx)
nhy <- nearest(rd$y, hy); dhy <- abs(rd$y - nhy)
on_v <- dvx <= dhy
line_coord <- ifelse(on_v, nvx, nhy)
pos        <- ifelse(on_v, rd$y, rd$x)
pos_orig   <- ifelse(on_v, y_min, x_min)
seg_idx <- pmin(pmax(floor((pos - pos_orig) / spacing) + 1L, 1L), n_seg_line)

key_rd  <- paste(ifelse(on_v, "v", "h"), round(line_coord, 1), seg_idx)
key_net <- paste(net$orientation, round(net$line_coord, 1), net$seg_idx)
rd$seg_id <- net$seg_id[match(key_rd, key_net)]
rd$class  <- net$class[rd$seg_id]
stopifnot(!anyNA(rd$class))

cls_lv <- c("freeway", "arterial", "collector", "local")

## ---------------------------------------------------------------------------
## 2. Pre-join diagnostics (the paper's audit, on the synthetic data)
## ---------------------------------------------------------------------------
cat("\n=== Diagnostics ===\n")
share_read <- 100 * prop.table(table(factor(rd$class, cls_lv)))
share_net  <- 100 * prop.table(table(factor(net$class, cls_lv)))  # equal-length segments
cat("Share of readings vs share of network length (%):\n")
print(round(rbind(readings = share_read, network = share_net), 1))

visits <- unique(rd[, c("seg_id", "date")])
rev_days <- tapply(visits$date, visits$seg_id, function(d) length(unique(d)))
seg_class <- net$class[as.integer(names(rev_days))]
med_rev <- tapply(rev_days, factor(seg_class, cls_lv), median)
cat("\nMedian revisit days per 100 m cell, by class (paper: 65/33/26/16):\n")
print(round(med_rev, 0))

## class means used to value the composition (computed from the demo data,
## as the paper values the blind-assigned mix at citywide class means)
rd$log_nox <- log(rd$NO + rd$NO2)
rd$log_pn  <- log(rd$PN1 + rd$PN2 + rd$PN3 + rd$PN4 + rd$PN5)
mu_nox <- tapply(rd$log_nox, factor(rd$class, cls_lv), mean)
mu_o3  <- tapply(rd$O3,      factor(rd$class, cls_lv), mean)
mu_pn  <- tapply(rd$log_pn,  factor(rd$class, cls_lv), mean)
cat("\nClass means: log(NOx), O3 (ppb), log(PN):\n")
print(round(rbind(log_nox = mu_nox, o3 = mu_o3, log_pn = mu_pn), 2))

## ---------------------------------------------------------------------------
## 3. Synthetic-outcome experiment: 500 non-freeway sites, blind vs constrained
## ---------------------------------------------------------------------------
n_sites <- 500L
nf <- net[net$class != "freeway", ]
seg_pick <- nf[sample(nrow(nf), n_sites, replace = TRUE), ]   # length-weighted (equal lengths)
u <- runif(n_sites)
site_x <- seg_pick$x0 + u * (seg_pick$x1 - seg_pick$x0)
site_y <- seg_pick$y0 + u * (seg_pick$y1 - seg_pick$y0)
site_class <- seg_pick$class

days_obs  <- sort(unique(rd$date))
hours_obs <- sort(unique(as.integer(substr(rd$Date_Time, 12, 13))))
site_day  <- sample(days_obs, n_sites, replace = TRUE)
site_hour <- sample(hours_obs, n_sites, replace = TRUE)
site_t <- as.numeric(as.POSIXct(paste0(site_day, " 00:00:00"), tz = "UTC")) +
  site_hour * 3600 + floor(runif(n_sites) * 3600)

## join metric: d = sqrt((ds/500 m)^2 + (dt/6 h)^2), caps 2000 m, +/-72 h
ord <- order(rd$t)
rt <- rd$t[ord]; rx <- rd$x[ord]; ry <- rd$y[ord]
rcl <- rd$class[ord]; rnox <- rd$log_nox[ord]; ro3 <- rd$O3[ord]; rpn <- rd$log_pn[ord]

match_one <- function(i) {
  lo <- findInterval(site_t[i] - 72 * 3600, rt) + 1L
  hi <- findInterval(site_t[i] + 72 * 3600, rt)
  if (hi < lo) return(c(NA, NA))
  idx <- lo:hi
  ds <- sqrt((rx[idx] - site_x[i])^2 + (ry[idx] - site_y[i])^2)
  keep <- ds <= 2000
  if (!any(keep)) return(c(NA, NA))
  idx <- idx[keep]
  d <- sqrt((ds[keep] / 500)^2 + ((abs(rt[idx] - site_t[i]) / 3600) / 6)^2)
  blind <- idx[which.min(d)]
  own <- rcl[idx] == site_class[i]
  constr <- if (any(own)) idx[own][which.min(d[own])] else NA
  c(blind, constr)
}
m <- vapply(seq_len(n_sites), match_one, numeric(2))
blind <- m[1, ]; constr <- m[2, ]

ok_b <- !is.na(blind); ok_p <- ok_b & !is.na(constr)
assigned_class <- rcl[blind[ok_b]]

## per-site quantities (composition gaps valued at demo class means, as in paper)
d_nox_comp <- mu_nox[assigned_class] - mu_nox[site_class[ok_b]]
d_o3_comp  <- mu_o3[assigned_class]  - mu_o3[site_class[ok_b]]
d_pn_comp  <- mu_pn[assigned_class]  - mu_pn[site_class[ok_b]]
d_nox_pair <- rnox[blind[ok_p]] - rnox[constr[ok_p]]
d_o3_pair  <- ro3[blind[ok_p]]  - ro3[constr[ok_p]]
d_pn_pair  <- rpn[blind[ok_p]]  - rpn[constr[ok_p]]

fw_share  <- 100 * mean(assigned_class == "freeway")
own_share <- 100 * mean(assigned_class == site_class[ok_b])

## ---------------------------------------------------------------------------
## 4. Bootstrap over sites (B = 199, percentile 95% CIs)
## ---------------------------------------------------------------------------
B <- 199L
nb <- sum(ok_b)
boot <- t(vapply(seq_len(B), function(b) {
  s <- sample.int(nb, nb, replace = TRUE)
  c(mean(d_nox_comp[s]), mean(d_o3_comp[s]), mean(d_pn_comp[s]))
}, numeric(3)))
np <- sum(ok_p)
bootp <- vapply(seq_len(B), function(b) {
  s <- sample.int(np, np, replace = TRUE)
  mean(d_nox_pair[s])
}, numeric(1))
ci <- function(v) quantile(v, c(0.025, 0.975))
ci_nox <- ci(boot[, 1]); ci_o3 <- ci(boot[, 2]); ci_pn <- ci(boot[, 3])
ci_pairnox <- ci(bootp)

## ---------------------------------------------------------------------------
## 5. Results summary
## ---------------------------------------------------------------------------
fmt_ci <- function(est, ci, d = 2)
  sprintf(paste0("%+.", d, "f  [%+.", d, "f, %+.", d, "f]"), est, ci[1], ci[2])

cat("\n=== Synthetic-outcome experiment (", n_sites, " non-freeway sites) ===\n", sep = "")
cat(sprintf("Blind-matched sites: %d (%.1f%%); paired (both matches): %d\n",
            nb, 100 * nb / n_sites, np))
cat(sprintf("Freeway-assigned share (blind): %.1f%%   Own-class share: %.1f%%\n",
            fw_share, own_share))
cat("\nComposition-implied gaps (blind mix vs own-class mix, valued at class means)\n")
cat("with B = 199 percentile bootstrap 95% CIs over sites:\n")
cat(sprintf("  log(NOx) inflation : %s log points\n", fmt_ci(mean(d_nox_comp), ci_nox)))
cat(sprintf("  O3 gap             : %s ppb\n",        fmt_ci(mean(d_o3_comp), ci_o3)))
cat(sprintf("  log(PN) gap        : %s log points\n", fmt_ci(mean(d_pn_comp), ci_pn)))
cat("\nPaired blind-vs-constrained gaps (same sites, same metric):\n")
cat(sprintf("  log(NOx): %s   O3: %+.2f ppb   log(PN): %+.2f\n",
            fmt_ci(mean(d_nox_pair), ci_pairnox), mean(d_o3_pair), mean(d_pn_pair)))

cat("\n=== Qualitative check against the paper's artifact ===\n")
chk <- c(
  "positive log(NOx) composition inflation" = mean(d_nox_comp) > 0,
  "negative (inverted) O3 composition gap"  = mean(d_o3_comp)  < 0,
  "small positive log(PN) composition gap"  = mean(d_pn_comp)  > 0 &&
                                              mean(d_pn_comp)  < mean(d_nox_comp),
  "revisit ordering local<collector<arterial<freeway" =
    all(diff(med_rev[c("local", "collector", "arterial", "freeway")]) > 0))
for (nm in names(chk)) cat(sprintf("  [%s] %s\n", ifelse(chk[nm], "PASS", "FAIL"), nm))
cat(ifelse(all(chk),
  "\nAll checks pass: the sampling-composition artifact reproduces qualitatively.\n",
  "\nWARNING: not all qualitative checks passed.\n"))
cat(sprintf("\nElapsed: %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
cat("(Synthetic demonstration only; magnitudes are not the paper's estimates.)\n")
