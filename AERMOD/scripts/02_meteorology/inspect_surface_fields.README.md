# inspect surface fields

Inspect surface meteorological fields for missing or invalid values.

## Inputs and outputs

Paths are specified in config/settings.json for configured stages. Population and statistics stages use TRIMAQS_DATA_ROOT or the command-line roots documented in their folder README. Read the settings block at the start of this file for its exact tables, column names and output paths. Input data are distributed separately.

## Use

Run the corresponding stage through run.py when it is listed in stages.json. Standalone population scripts are run as saved Python files; statistics scripts accept --help. This file may be a helper imported by another script, rather than a separate stage.

## Assumptions and validation

Retain the units and definitions in the dataset READMEs. Use the final selected inventories and common grid. Syntax checks do not establish scientific agreement; verify calendars, mass balance, grid alignment and category sums before using new results. Preserve existing results when changing settings.

## Configured files

- `inspect_surface_fields.boundary_input_dir_5`: `inputs/meteorology/era5/Raw files/2024 data for time correction`
- `inspect_surface_fields.cloud_input_dir_4`: `inputs/meteorology/era5/Raw files/ERA5 surface data - raw/Single variables`
- `inspect_surface_fields.grid_csv_file_1`: `inputs/supporting/grids/meteorological_grid_10km.csv`
- `inspect_surface_fields.land_input_dir_3`: `inputs/meteorology/era5/Raw files/ERA5 surface data - raw/zipped data`
- `inspect_surface_fields.surface_base_input_dir_2`: `inputs/meteorology/era5/Raw files/ERA5 surface data - raw/ERA5`
- `inspect_surface_fields.surface_output_dir_6`: `inputs/meteorology/era5/ERA5 surface yearly/rebuilt_utc`
- `inspect_surface_fields.upper_input_dir_8`: `inputs/meteorology/era5/Raw files/ERA5 upper air data/pressure_level_data`
- `inspect_surface_fields.upper_output_dir_7`: `inputs/meteorology/era5/ERA5 upper air yearly/rebuilt_utc`
- `surface_inspection_report`: `working/meteorology/surface_inspection.json`
