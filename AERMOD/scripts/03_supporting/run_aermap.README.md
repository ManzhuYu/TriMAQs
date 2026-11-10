# run aermap

Prepare terrain and receptor inputs with AERMAP.

## Inputs and outputs

Paths are specified in config/settings.json for configured stages. Population and statistics stages use TRIMAQS_DATA_ROOT or the command-line roots documented in their folder README. Read the settings block at the start of this file for its exact tables, column names and output paths. Input data are distributed separately.

## Use

Run the corresponding stage through run.py when it is listed in stages.json. Standalone population scripts are run as saved Python files; statistics scripts accept --help. This file may be a helper imported by another script, rather than a separate stage.

## Assumptions and validation

Retain the units and definitions in the dataset READMEs. Use the final selected inventories and common grid. Syntax checks do not establish scientific agreement; verify calendars, mass balance, grid alignment and category sums before using new results. Preserve existing results when changing settings.

## Configured files

- `run_aermap.coordinate_templates`: `inputs/supporting/receptors`
- `run_aermap.dem_folder`: `inputs/supporting/dem`
- `run_aermap.executable`: `executables/aermap/aermap`
- `run_aermap.grid_file`: `inputs/supporting/grids/meteorological_grid_10km.csv`
- `run_aermap.output_folder`: `working/aermap`
- `run_aermap.temporary_dem_folder`: `working/dem_tiles`
