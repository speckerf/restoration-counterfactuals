# Restoration counterfactuals

Two sequential notebooks select comparable donor pixels for one restoration site, then estimate separate ridge synthetic controls for short vegetation height (SVH) and aboveground biomass (AGB).

1. [Select controls](notebooks/01_select_controls.ipynb): authenticate Earth Engine, read the boundary and `config.yaml`, sample treatment pixels, exclude spillover and other projects, inspect matching masks, save locations.
2. [Estimate impacts](notebooks/02_estimate_impacts.ipynb): extract annual boundary means and donor-pixel values, fit each indicator, save time series, effects, weights, diagnostics and PDF/PNG figures in `outputs/`.

Install locally with `uv sync --locked`, select `.venv/bin/python` as your notebook kernel, and run notebook 1 then 2. Alternatively use `python -m pip install -e .` with Python 3.11+. Functions live in `src/restoration_counterfactuals`; notebooks handle sequencing and saving.

For [Google Colab](https://colab.research.google.com/), upload each notebook using **File → Upload notebook** and upload/unzip the repository into `/content/restoration-counterfactuals`. Once this repository has a public Git URL, set `REPO_URL` in the installation cell instead. No Git remote is configured yet, so repository-specific Colab links cannot be supplied. Each notebook installs from `pyproject.toml`. Download `outputs/matching_outputs.zip` after notebook 1 and extract it into the repository root in notebook 2's runtime. Set `EE_PROJECT` and authenticate in each runtime.

## Examples and exclusion coverage

| GeoJSON | Intervention year |
| --- | --- |
| `examples/site.geojson` / `finca_rinconada.geojson` (default) | 2022 |
| `examples/campeche_mexico.geojson` | 2021 |
| `examples/bom_futuro.geojson` | 2017 |
| `examples/reforest_now_au.geojson` | 2020 |

Coordinates and years were supplied by the user. The supplied “finca rotoanda” year is assigned to the `finca_rinconada` polygon. Change `SITE_PATH` in notebook 1 to choose another site; notebook 2 reads the exact site saved by notebook 1.

`examples/exclusion_boundaries.geojson` is an **empty FeatureCollection at user request**. No neighboring restoration projects are excluded for this run. The donor-search function retains buffered exclusion support: add boundaries to this file when available. No full Restor database is required. Both notebooks default to the supplied Earth Engine project `ee-speckerfelix`.

## Data access and definitions

Use an Earth Engine account and a Cloud project registered for Earth Engine, with permission to use its API. The loaders reference these hosted datasets; access and coverage must be available to the authenticated account:

- AlphaEarth: `GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL`, intervention year minus one for interventions from 2018 ([catalog](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL)).
- Earlier interventions: GLAD `projects/glad/GLCLU2020/v2/LCLUC_{year}`, latest five-year snapshot strictly before intervention, remapped to broad classes; class 255 excluded. Interventions must be after 2000.
- Human modification: `projects/sat-io/open-datasets/GHM/HM_1990_2020_OVERALL_300M`, interpolation at intervention year minus one, capped at 2020 as in the source.
- Disturbance: `UMD/hansen/global_forest_change_2025_v1_13`, pre-intervention loss years. Zero loss and later loss are normalized to 9999, correcting the source loader/matcher encoding mismatch. Genuine missing coverage follows the original unmask behavior.
- SVH: `projects/global-pasture-watch/assets/gsvh-30m/v1/short-veg-height_m`, multiplied by 0.1, metres.
- AGB: `projects/sat-io/open-datasets/CTREES-GLOBAL-AGB-100M`, band `agb`, multiplied by 0.1, Mg/ha; negative values masked.

Annual years come from source image IDs. Treated outcomes are raster means directly inside the project boundary, not averages of treatment samples. Both polygon reduction and donor sampling use each outcome's native grid ([Earth Engine reduction reference](https://developers.google.com/earth-engine/apidocs/ee-image-reduceregions)). Missing annual observations remain missing. Distinct matching-grid donor coordinates may still fall in the same coarser AGB pixel.

## Fit and validation

The Python solver ports `R/utils_scm_fit.R`: nonnegative weights summing to one, raw univariate outcomes, pre-intervention observations only, and `D = 2(XᵀX + λI) + 1e-8 I`. It uses the same quadprog algorithm, truncates weights below 1e-10 and renormalizes. Of 30 log-spaced penalties from 1e-6 to 1e3, select the largest within 2% of minimum pre-fit RMSPE. Exact zero error is handled without dividing by zero.

At least three complete pre-intervention years and one donor are required. Missing years are excluded from fitting, retained as gaps in trajectories, and counted in diagnostics. Duplicate sampled coordinates are collapsed so repetition does not alter ridge regularization. No multivariate fitting, other estimators, or uncertainty intervals are included. Pre-fit agreement alone does not establish causality.

Run `uv run pytest`. The committed numerical fixture was evaluated with the original R `solve_scm_ridge` function and its penalty-selection rule; tests compare weights, penalty and RMSPE, and check missingness and pre-only fitting. The notebook computation cells also ran successfully against live Earth Engine using `ee-speckerfelix`: Finca Rinconada, 80 unique donors, 22 complete pre-intervention years per indicator. Both live fits were compared against the original R solver with `tests/compare_r_fit.R`; weights agreed within 1e-10 and selected penalties matched. The Bom Futuro land-cover lookup and undisturbed matching branch were checked live. Generated results are in the gitignored `outputs/` directory. Fresh Colab execution remains unverified because the available browser session is signed out; installation/authentication cells were not exercised there.
