# Input data

All locations below are relative to `data_root` in the configuration.

| Folder | Contents |
|---|---|
| `inputs/emissions/drilling_completion` | Daily drilling and completion emissions, tonnes per day |
| `inputs/emissions/production` | Combined leakage and heater emissions, tonnes per month |
| `inputs/emissions/compressors` | Compressor emissions, tonnes per month |
| `inputs/emissions/other_tier2` | Other Tier 2 facility emissions, tonnes per month |
| `inputs/parameters` | Source locations, elevations, stack parameters and case selections |
| `inputs/supporting` | Grids, study boundary, receptors, receptor mapping and raster mask |
| `inputs/meteorology` | ERA5, AERSURFACE and AERMET inputs/products needed for model runs |
| `inputs/ml_predictions` | Daily ML predictions used for source contributions |
| `inputs/source_data` | Optional upstream inputs for inventory construction |

Inventory names follow `compressor_emissions_ch4_2007.csv` or `other_tier2_emissions_nox_2007.csv`, within category/pollutant folders. Every inventory uses `source_id`. Daily fields are ISO dates (`YYYY-MM-DD`); monthly fields are `jan_tonnes` through `dec_tonnes`. NO2 modeling uses the NOX inventory; it is not an additional emission inventory to sum with NOX.

The parameter folder contains four location tables: `drilling_completion_data.csv`, `production_data.csv`, `compressor_data_with_elevation.csv` and `other_tier2_data_with_elevation.csv`. `stack_parameters.csv` gives point-source stack parameters, area-source dimensions and allocation weights. The source_type field identifies the applicable fields. `model_cases.csv` connects each source/year/pollutant to its location, releases and meteorological/receptor grids. Identifiers are text. `location_id` distinguishes locations where a source has different case settings.

Use `data_dictionary.csv` and the folder READMEs supplied with the inputs for field definitions. Heights, elevations and dimensions use metres; temperatures use kelvin; exit velocities use metres per second. Monthly masses are divided by calendar days before hourly allocation. Production already includes leakage and heater emissions.

The final tables are sufficient for the dispersion stages. Rebuilding inventories from their original sources requires the separate activity records, factors and worksheets. Large meteorological files and model executables may be supplied separately.
