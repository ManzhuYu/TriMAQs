# run other tier2

Run AERMOD for other UOGD facilities.

## Inputs and outputs

Paths are specified in config/settings.json for configured stages. Population and statistics stages use TRIMAQS_DATA_ROOT or the command-line roots documented in their folder README. Read the settings block at the start of this file for its exact tables, column names and output paths. Input data are distributed separately.

## Use

Run the corresponding stage through run.py when it is listed in stages.json. Standalone population scripts are run as saved Python files; statistics scripts accept --help. This file may be a helper imported by another script, rather than a separate stage.

## Assumptions and validation

Retain the units and definitions in the dataset READMEs. Use the final selected inventories and common grid. Syntax checks do not establish scientific agreement; verify calendars, mass balance, grid alignment and category sums before using new results. Preserve existing results when changing settings.

## Configured files

- `run_other_tier2.aermap_gridwise_dir_4`: `inputs/supporting/receptors`
- `run_other_tier2.aermod_executable_1`: `executables/aermod/aermod`
- `run_other_tier2.grid_csv_file_3`: `inputs/supporting/grids/meteorological_grid_10km.csv`
- `run_other_tier2.inventory_template_7`: `inputs/emissions/other_tier2/{pollutant}/other_tier2_emissions_{pollutant_slug}_{year}.csv`
- `run_other_tier2.log_directory_6`: `outputs/model_runs/logs/other_tier2`
- `run_other_tier2.meteorological_data_dir_template_2`: `inputs/meteorology/aermet/{year}/Grid_{grid_ID}`
- `run_other_tier2.output_base_dir_template_5`: `outputs/model_runs/other_tier2/{year}/{pollutant}/{facility_id}`
