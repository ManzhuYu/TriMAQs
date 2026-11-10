# prepare era5

Prepare annual surface and upper-air meteorological inputs from raw ERA5 files.

## Inputs and outputs

Paths are specified in config/settings.json for configured stages. Population and statistics stages use TRIMAQS_DATA_ROOT or the command-line roots documented in their folder README. Read the settings block at the start of this file for its exact tables, column names and output paths. Input data are distributed separately.

## Use

Run the corresponding stage through run.py when it is listed in stages.json. Standalone population scripts are run as saved Python files; statistics scripts accept --help. This file may be a helper imported by another script, rather than a separate stage.

## Assumptions and validation

Retain the units and definitions in the dataset READMEs. Use the final selected inventories and common grid. Syntax checks do not establish scientific agreement; verify calendars, mass balance, grid alignment and category sums before using new results. Preserve existing results when changing settings.

## Configured files

- `prepare_era5.boundary_input_dir_5`: `inputs/meteorology/era5/Raw files/2024 data for time correction`
- `prepare_era5.cloud_input_dir_4`: `inputs/meteorology/era5/Raw files/ERA5 surface data - raw/Single variables`
- `prepare_era5.grid_csv_file_1`: `inputs/supporting/grids/meteorological_grid_10km.csv`
- `prepare_era5.land_input_dir_3`: `inputs/meteorology/era5/Raw files/ERA5 surface data - raw/zipped data`
- `prepare_era5.surface_base_input_dir_2`: `inputs/meteorology/era5/Raw files/ERA5 surface data - raw/ERA5`
- `prepare_era5.surface_output_dir_6`: `inputs/meteorology/era5/ERA5 surface yearly/rebuilt_utc`
- `prepare_era5.upper_input_dir_8`: `inputs/meteorology/era5/Raw files/ERA5 upper air data/pressure_level_data`
- `prepare_era5.upper_output_dir_7`: `inputs/meteorology/era5/ERA5 upper air yearly/rebuilt_utc`
