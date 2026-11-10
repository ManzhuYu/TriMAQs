# allocate development emissions

Allocate drilling and completion emissions to daily activity intervals.

## Inputs and outputs

Paths are specified in config/settings.json for configured stages. Population and statistics stages use TRIMAQS_DATA_ROOT or the command-line roots documented in their folder README. Read the settings block at the start of this file for its exact tables, column names and output paths. Input data are distributed separately.

## Use

Run the corresponding stage through run.py when it is listed in stages.json. Standalone population scripts are run as saved Python files; statistics scripts accept --help. This file may be a helper imported by another script, rather than a separate stage.

## Assumptions and validation

Retain the units and definitions in the dataset READMEs. Use the final selected inventories and common grid. Syntax checks do not establish scientific agreement; verify calendars, mass balance, grid alignment and category sums before using new results. Preserve existing results when changing settings.

## Configured files

- `allocate_development_emissions.completion_workbook`: `inputs/source_data/activity/completion_emissions.xlsx`
- `allocate_development_emissions.drilling_workbook`: `inputs/source_data/activity/drilling_emissions.xlsx`
- `allocate_development_emissions.pad_mapping`: `inputs/parameters/well_to_pad.csv`
- `allocate_development_emissions.pad_output`: `working/inventory/development_by_pad`
- `allocate_development_emissions.transport_workbook`: `inputs/source_data/activity/transport_emissions.xlsx`
- `allocate_development_emissions.well_output`: `working/inventory/development_by_well`
