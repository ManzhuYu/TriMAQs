# run aermet

Run both AERMET stages for selected years and meteorological grids.

## Inputs and outputs

Paths are specified in config/settings.json for configured stages. Population and statistics stages use TRIMAQS_DATA_ROOT or the command-line roots documented in their folder README. Read the settings block at the start of this file for its exact tables, column names and output paths. Input data are distributed separately.

## Use

Run the corresponding stage through run.py when it is listed in stages.json. Standalone population scripts are run as saved Python files; statistics scripts accept --help. This file may be a helper imported by another script, rather than a separate stage.

## Assumptions and validation

Retain the units and definitions in the dataset READMEs. Use the final selected inventories and common grid. Syntax checks do not establish scientific agreement; verify calendars, mass balance, grid alignment and category sums before using new results. Preserve existing results when changing settings.

## Configured files

- `run_aermet.aermap_output_dir_5`: `inputs/supporting/receptors`
- `run_aermet.aermet_executable_6`: `executables/aermod/aermet`
- `run_aermet.aersurface_output_dir_4`: `inputs/meteorology/aersurface`
- `run_aermet.grid_csv_file_1`: `inputs/supporting/grids/meteorological_grid_10km.csv`
- `run_aermet.output_dir_7`: `inputs/meteorology/aermet`
- `run_aermet.surface_nc_file_folder_2`: `inputs/meteorology/era5/ERA5 surface yearly/surface_complete`
- `run_aermet.upperair_nc_file_folder_3`: `inputs/meteorology/era5/ERA5 upper air yearly/rebuilt_utc`
