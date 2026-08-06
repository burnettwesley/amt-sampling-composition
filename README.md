# Code for: Sampling-composition bias in urban mobile air-quality monitoring

Replication code for:

> Burnett, J. Wesley and Tamara L. Sheldon (2026). "Sampling-composition bias in
> urban mobile air-quality monitoring: diagnostics and a road-class-aware
> correction for proximity joins to located outcomes." Manuscript prepared for
> *Atmospheric Measurement Techniques*.

All analysis is in R. Scripts are numbered as in the authors' working
repository and should be run from the repository root, in the order below.
No monitoring data are redistributed here; see "Data access" for how to obtain
the inputs.

## Demo (no data agreement required)

A fully synthetic demonstration reproduces the paper's mechanism end-to-end,
offline, in seconds, with base R only:

```
Rscript scripts/00_generate_demo_data.R   # build synthetic network + ~200k readings
Rscript scripts/99_demo_run.R             # diagnostics + join experiment + bootstrap
```

The first script builds a stylized 30x30 street grid with a freeway spine and
simulates 120 collection days of two-car routes calibrated to the paper's
published revisit-day medians (65/33/26/16) and class-level concentration
means, writing `demo_data/demo_readings.csv.gz` (~7 MB). The second classifies
every reading by road class, runs the pre-join diagnostics, matches 500
synthetic non-freeway sites blind versus class-constrained under the paper's
join metric (2 km, +/-72 h), and bootstraps the gaps (B = 199). The artifact
reproduces qualitatively: blind joins inflate log(NOx) at non-freeway sites,
understate the titrated pollutant O3 (negative gap), and carry a small
positive particle-number gap tracking that pollutant's shallow gradient.

**Disclaimer (loud, on purpose): the demo data are SYNTHETIC. They contain no
Aclima records and were generated exclusively from summary statistics
published in the manuscript; the generator reads no Aclima file. Structural
properties only — not for scientific inference.** See `demo_data/README.md`
for the parameterization and compliance basis.

## Script run order and outputs

| Order | Script                            | What it does and produces |
|------:|-----------------------------------|---------------------------|
|     1 | `33_reading_road_class.R`         | Parses and classifies the Oakland OpenStreetMap network (freeway / arterial / collector / local) and snaps every Oakland reading to its nearest segment. Produces `data/processed/osm_roads.gpkg` and `data/processed/reading_road_class.csv.gz`. Steps 3 and 4 classify companion-project records and are not required for this paper. |
|     2 | `44_sf_osm.R`                     | Downloads the San Francisco OSM network via the Overpass API, applies the identical class mapping, and classifies every SF reading. Produces `data/processed/sf_osm_roads.gpkg`, `data/processed/sf_reading_road_class.csv.gz`, and `data/processed/sf_network_km.csv`. |
|     3 | `42_cross_car_agreement.R`        | Between-vehicle collocation check on Oakland dual-coverage cell-hours (n = 1,198). Produces `output/tables/cross_car_agreement.txt`. |
|     4 | `54_amt_coverage_map.R`           | Two-panel revisit-frequency coverage map (paper Fig. 1). Produces `output/figures/fig_amt_coverage.png`. |
|     5 | `45_sf_siting.R`                  | Composition diagnostics (readings vs. network road-km, revisit frequency, class means) and the archived synthetic-site matching experiment for both cities. Site hours in this archived version are drawn from a companion-project file; script 56 is the outcome-free rerun reported in the paper. Produces `output/tables/sf_siting.txt`. |
|     6 | `55_amt_sensitivity.R`            | Parameter sweep over the join caps, cell-size robustness, and revisit distributions (paper Fig. 2 and Table 1). Produces `output/tables/amt_sensitivity.txt` and `output/figures/fig_amt_revisit_dist.png`. |
|     7 | `56_synthetic_uniform_hours.R`    | Outcome-free synthetic-site experiment with uniform collection-hour draws; the paper's baseline numbers. Produces `output/tables/amt_synthetic_uniform.txt`. |
|     8 | `64_o3_sign_reversal.R`           | Ozone sign-reversal replication (San Francisco; identical sites and matches, O3 read off instead of NOx). Produces `output/tables/o3_sign_reversal.txt`. |
|     9 | `65_ufp_replication.R`            | Particle-number (UFP) replication (San Francisco; identical sites and matches). Produces `output/tables/ufp_replication.txt`. |
|    10 | `67_amt_bootstrap.R`              | Nonparametric bootstrap over synthetic sites (B = 999, percentile 95% intervals) plus the cross-car correlation interval; reproduces the scripts 56/64/65 baselines before any interval is computed. Produces `output/tables/amt_bootstrap.txt`. |
|    11 | `68_distance_decay.R`             | Distance-decay of NOx away from the freeway network, non-freeway readings only, both cities. Produces `output/tables/distance_decay.txt`. |

## Expected directory layout

Scripts expect the following relative to the repository root:

```
data/raw/aclima/Oakland_201505_201605_GoogleAclimaAQ.txt       (Oakland release)
data/raw/aclima/California_201605_201709_GoogleAclimaAQ.txt    (California release)
data/raw/osm/                                                  (cached Overpass JSON)
data/processed/                                                (script outputs, created)
output/tables/  output/figures/                                (script outputs, created)
```

The California release path can instead be supplied through the `AQ_CA_FILE`
environment variable.

## Data access

- **Aclima/EDF mobile monitoring data** (Google Street View--Aclima campaigns,
  Oakland 2015--2016 and California/San Francisco 2016--2017): available through
  the Environmental Defense Fund air-quality data portal; registration is
  required. No monitoring data are redistributed in this repository.
- **OpenStreetMap road networks**: retrieved via the Overpass API
  (`44_sf_osm.R` performs the San Francisco download; the Oakland network JSON
  is cached the same way). OSM data are (c) OpenStreetMap contributors, ODbL.

## License

MIT; see `LICENSE`. Please cite the paper (see `CITATION.cff`) if you use this
code.
