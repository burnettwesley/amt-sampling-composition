# DESTINATION: AMT paper (see docs/SPLIT_ALLOCATION.md)
# 42_cross_car_agreement.R
# ========================
# Oakland-specific between-vehicle comparability check (Wes's suggestion,
# 2026-07-19). The two Street View cars rarely sampled the same place at the
# same time (divided polygon assignments), but in the cell-hours where both
# did, their readings provide an in-sample collocation test -- the paper
# currently borrows all between-vehicle evidence from Solomon et al. (2020),
# whose campaigns ran in other cities.
#
# Design: ~100m grid cells x clock hour. For cell-hours where both cars
# recorded valid NOx, compare per-car mean log(NOx): correlation, regression
# slope, and Solomon-style fractional absolute mean difference (FAMD =
# |A-B| / mean(A,B), in levels). Benchmark against a WITHIN-car split-half
# FAMD (odd vs even readings of one car in the same cell-hour), which bounds
# what agreement is achievable given within-hour temporal variability alone.
# Also: road-class composition of the dual cell-hours, and BC availability
# by car (Solomon reports BC rode on one vehicle at a time in their campaign).
#
# OUTPUT: output/tables/cross_car_agreement.txt
# =============================================================================

suppressPackageStartupMessages({library(tidyverse); library(lubridate)})

a <- read_csv("data/raw/aclima/Oakland_201505_201605_GoogleAclimaAQ.txt",
              col_select = c(Date_Time, Car_Identifier, Latitude, Longitude,
                             NO2, NO, BC),
              show_col_types = FALSE) %>%
  mutate(ts  = ymd_hms(Date_Time, tz = "UTC"),
         NO2 = suppressWarnings(as.numeric(NO2)),
         NO  = suppressWarnings(as.numeric(NO)),
         BC  = suppressWarnings(as.numeric(BC))) %>%
  filter(!is.na(Latitude), !is.na(Longitude), !is.na(ts))

rc <- read_csv("data/processed/reading_road_class.csv.gz", show_col_types = FALSE)
stopifnot(nrow(rc) == nrow(a))
a$road_class <- rc$road_class

a <- a %>%
  mutate(cell = paste(floor(Longitude / 0.0012), floor(Latitude / 0.0009)),
         hr   = floor_date(ts, "hour"),
         nox  = NO2 + NO) %>%
  filter(!is.na(nox), nox > 0)

MIN_READS <- 5   # per car per cell-hour, for stable means

# -- Cross-car cell-hour panel ------------------------------------------------
cc <- a %>%
  group_by(cell, hr, Car_Identifier) %>%
  summarise(n = n(), mean_nox = mean(nox),
            class = names(sort(table(road_class), decreasing = TRUE))[1],
            .groups = "drop") %>%
  filter(n >= MIN_READS) %>%
  pivot_wider(id_cols = c(cell, hr), names_from = Car_Identifier,
              values_from = c(mean_nox, n, class)) %>%
  filter(!is.na(mean_nox_Car_A), !is.na(mean_nox_Car_B)) %>%
  mutate(lA = log(mean_nox_Car_A), lB = log(mean_nox_Car_B),
         famd = abs(mean_nox_Car_A - mean_nox_Car_B) /
                ((mean_nox_Car_A + mean_nox_Car_B) / 2))

# -- Within-car split-half benchmark (same cell-hours' populations) -----------
sh <- a %>%
  group_by(cell, hr, Car_Identifier) %>%
  filter(n() >= 2 * MIN_READS) %>%
  mutate(half = row_number() %% 2) %>%
  group_by(cell, hr, Car_Identifier, half) %>%
  summarise(m = mean(nox), .groups = "drop") %>%
  pivot_wider(names_from = half, values_from = m, names_prefix = "h") %>%
  filter(!is.na(h0), !is.na(h1)) %>%
  mutate(famd = abs(h0 - h1) / ((h0 + h1) / 2))

dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
sink("output/tables/cross_car_agreement.txt")
cat("Cross-car agreement check -- scripts/42 --", format(Sys.Date()), "\n")
cat(sprintf("Cell = ~100m grid; hour = clock hour; MIN %d readings/car/cell-hour.\n\n",
            MIN_READS))

cat("== BC availability by car (share of readings with non-missing BC) ==\n")
bc <- a %>% group_by(Car_Identifier) %>%
  summarise(n = n(), bc_share = 100 * mean(!is.na(BC)), .groups = "drop")
for (i in seq_len(nrow(bc)))
  cat(sprintf("  %s: %.1f%% of %s readings\n", bc$Car_Identifier[i],
              bc$bc_share[i], format(bc$n[i], big.mark = ",")))
cat("\n")

cat("== Dual-coverage cell-hours (both cars, >=5 valid NOx readings each) ==\n")
cat(sprintf("n = %s cell-hours\n", format(nrow(cc), big.mark = ",")))
comp <- cc %>% count(class_Car_A) %>% mutate(pct = 100 * n / sum(n)) %>%
  arrange(desc(pct))
cat("road-class composition (by Car_A readings' modal class):\n")
for (i in seq_len(nrow(comp)))
  cat(sprintf("  %-9s %5.1f%%\n", comp$class_Car_A[i], comp$pct[i]))
cat("\n")

cat("== Cross-car agreement, mean log(NOx) per cell-hour ==\n")
cat(sprintf("correlation: %.3f\n", cor(cc$lA, cc$lB)))
fit <- lm(lB ~ lA, data = cc)
cat(sprintf("regression of B on A: slope %.3f (SE %.3f), R2 %.3f\n",
            coef(fit)[2], summary(fit)$coefficients[2, 2],
            summary(fit)$r.squared))
cat(sprintf("FAMD (levels): median %.1f%%, mean %.1f%%, share < 20%%: %.1f%%\n\n",
            100 * median(cc$famd), 100 * mean(cc$famd),
            100 * mean(cc$famd < 0.20)))

cat("== Within-car split-half benchmark (same-cell-hour halves, one car) ==\n")
cat(sprintf("n = %s car-cell-hours\n", format(nrow(sh), big.mark = ",")))
cat(sprintf("FAMD (levels): median %.1f%%, mean %.1f%%, share < 20%%: %.1f%%\n",
            100 * median(sh$famd), 100 * mean(sh$famd),
            100 * mean(sh$famd < 0.20)))
cat("\nReading: if cross-car FAMD is comparable to the within-car split-half\n")
cat("FAMD, between-vehicle differences are no larger than the within-hour\n")
cat("sampling variability a single car exhibits at the same location.\n")

# -- Same-car CROSS-PASS benchmark (the fair comparison) ----------------------
# Split-half pairs readings seconds apart; cross-car pairs can be up to an
# hour apart. The fair benchmark is the same car making two distinct passes
# through the same cell within the same hour (pass break = >120 s gap).
passes <- a %>%
  arrange(Car_Identifier, cell, hr, ts) %>%
  group_by(Car_Identifier, cell, hr) %>%
  mutate(gap = as.numeric(ts - lag(ts), units = "secs"),
         pass = cumsum(is.na(gap) | gap > 120)) %>%
  group_by(Car_Identifier, cell, hr, pass) %>%
  summarise(n = n(), m = mean(nox), t_mid = mean(ts), .groups = "drop") %>%
  filter(n >= MIN_READS) %>%
  group_by(Car_Identifier, cell, hr) %>%
  filter(n() >= 2) %>%
  arrange(t_mid) %>%
  summarise(m1 = m[1], m2 = m[2],
            sep_min = as.numeric(t_mid[2] - t_mid[1], units = "mins"),
            .groups = "drop") %>%
  mutate(famd = abs(m1 - m2) / ((m1 + m2) / 2))

# time separation of the cross-car pairs, for comparability
cc_sep <- a %>%
  group_by(cell, hr, Car_Identifier) %>%
  filter(n() >= MIN_READS) %>%
  summarise(t_mid = mean(ts), .groups = "drop") %>%
  pivot_wider(id_cols = c(cell, hr), names_from = Car_Identifier,
              values_from = t_mid) %>%
  filter(!is.na(Car_A), !is.na(Car_B)) %>%
  mutate(sep_min = abs(as.numeric(Car_A - Car_B, units = "mins")))

cat("\n== Same-car cross-pass benchmark (two passes, same cell-hour) ==\n")
cat(sprintf("n = %s car-cell-hours; median pass separation %.1f min\n",
            format(nrow(passes), big.mark = ","), median(passes$sep_min)))
cat(sprintf("FAMD (levels): median %.1f%%, mean %.1f%%, share < 20%%: %.1f%%\n",
            100 * median(passes$famd), 100 * mean(passes$famd),
            100 * mean(passes$famd < 0.20)))
cat(sprintf("correlation of log pass means: %.3f\n",
            cor(log(passes$m1), log(passes$m2))))
cat(sprintf("\ncross-car pairs' median time separation: %.1f min\n",
            median(cc_sep$sep_min)))
cat("\nReading: the same-car cross-pass FAMD is the fair benchmark for the\n")
cat("cross-car FAMD; if the two are similar at similar time separations,\n")
cat("the cross-car disagreement reflects within-hour pollution variability\n")
cat("rather than between-vehicle instrument differences.\n")
sink()
cat(readLines("output/tables/cross_car_agreement.txt"), sep = "\n")
