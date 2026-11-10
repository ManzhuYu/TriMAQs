# TriMAQs emission and concentration workflow

Code accompanying “Trends in Exposure to Air Pollutants from Marcellus and Utica Shale Development (TriMAQs)” by Gernand et al., covering 2004–2023.

## Setup

Install the packages in requirements.txt in a Python environment with the geospatial system libraries. AERMOD and AERMET version 24142 and the terrain/surface preprocessors are separate executables. Copy config/settings.example.json to config/settings.json and set data_root to your downloaded TriMAQs_data folder. All configured paths are relative to that root. Alternatively set TRIMAQS_DATA_ROOT. No HPC account or personal directory is required.

## Run order

1. scripts/01_inventory contains the activity, production and compressor calculations. Use the archived, selected model-ready inventories to reproduce the final study: do not add their raw components a second time. The generic facility interpolation script is an upstream utility, not a replacement for the final reviewed facility selection.
2. scripts/02_meteorology prepares ERA5 and AERMET inputs. Full meteorology is obtained separately if it is not in the data download.
3. scripts/03_supporting prepares terrain and the receptor mapping.
4. scripts/04_dispersion runs the four UOGD categories with the shared source-parameter reader. FASTALL is disabled. Final inventories, model_cases.csv, source locations and stack parameters are required.
5. scripts/05_postprocessing compiles receptor concentrations and writes concentration and source-contribution NetCDF files.
6. scripts/06_population rebuilds total and demographic populations from the supplied Census/NHGIS source archive if needed. The shared annual CSVs can be used directly.
7. scripts/07_statistics calculates current concentration and exposure summaries. See that folder’s README for the commands.
8. notebooks/Report_figures.ipynb creates figures, one figure per headed section. Set TRIMAQS_DATA_ROOT before starting Jupyter, or set DATA_ROOT and RESULTS_ROOT in its settings cell.

Run configured stages with `python run.py STAGE --config config/settings.json`. Use `--check` first to resolve configuration and check syntax. This check does not validate scientific results. Each Python file has a neighboring README with inputs, outputs, use and limitations.

## Data and interpretation

The public data folder contains inputs, outputs and supporting spatial data. Outputs/concentrations contains daily and annual fields; outputs/source_contributions contains daily ratios and two different annual summaries. NO2 and SO2 concentrations are ppb; primary PM2.5 and VOC are µg/m³. Tier_1 plus Tier_2 equals UOGD: do not sum all three source labels. NOX is retained as an emission-inventory pollutant; CH4 is not a report concentration output.

Source/ML ratios are descriptive comparisons between independent model estimates. Missing or nonpositive ML denominators stay missing; values above one are retained. The mean of daily ratios is not the ratio of annual means. Neither is a closed source budget.

Population files contain nonnegative integer estimates allocated by tract–grid overlap. ACS periods use midpoint years; 2004–2007 are interpolated and 2022–2023 extrapolated with nonnegative constraints. Low income means an income/poverty ratio below 2. The combined non-Hispanic AIAN/Asian/Pacific group is not Asian alone.

CNX comparisons, diagnostic reruns, machine-specific batch files and private checkpoints are not part of this repository. Existing uncertainty/sensitivity experiments are not presented as recomputed under the final configuration.

## Reproducibility status

This is a submission draft. The final HPC outputs have been summarized separately; the relocated package has syntax and focused synthetic checks. A complete execution of every stage from raw inputs through the portable release has not been performed. Upstream third-party datasets and executable availability must be checked before a new full run. The report and validation sections are still under review. No public data license is implied for third-party source material.

## Compact public data release

The public archive retains UOGD only in concentration and source-ratio NetCDFs (source dimension length one). Tier-specific outputs, embedded ML concentration fields, ML prediction files, receptor files, statistics and meteorology remain outside this download. Annual population CSVs and optional annual concentration CSVs are included. Reproducing all model stages requires the external inputs described above. Figure sections needing private tier outputs or statistics cannot run from the compact archive alone; the annual UOGD map reader accepts the one-source layout.
