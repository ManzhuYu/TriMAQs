# Emission inventories

These scripts support inventory preparation. To reproduce the final concentration dataset, start with the final inventories supplied in the data package. Upstream worksheets, activity data and emission factors are needed only when rebuilding the inventories.

| Script | Purpose |
|---|---|
| `estimate_activity_emissions.py` | General activity × emission-factor calculation with explicit input units |
| `allocate_development_emissions.py` | Allocate drilling and completion mass to daily activity periods |
| `aggregate_production_to_pads.py` | Combine well-level production emissions at the pad level |
| `build_production_inventory.py` | Combine leakage and heater inputs into monthly pad inventories |
| `estimate_compressor_emissions.py` | Estimate compressor emissions using production and factors |
| `interpolate_facility_emissions.py` | Interpolate facility emissions across reporting years |

Most stages use the shared configuration, for example:

```bash
python run.py build_production_inventory --config config/settings.json
```

The general activity calculator has its own command-line interface; use `python scripts/01_inventory/estimate_activity_emissions.py --help`. It is a unit-explicit calculation tool, not a replacement for the study's underlying parameter worksheets.

Drilling/completion model inputs are daily mass tables. Combined production inputs are monthly mass tables. NOx inventory values feed the separate NOX and NO2 model runs; those outputs must not be added together. The combined production inventory includes leakage and heating, with no production SO2 estimate. Check monthly/daily definitions and mass units in the source inventory dictionary before rebuilding compressor or facility inventories. These preparation scripts retain their source-table column expectations; not every upstream table is included in the final model-ready input package.

Outputs should go to a new preparation location. Review mass balances, well-to-pad matching and units before replacing any final inventory. File paths and available year/pollutant selections are listed in `config/settings.example.json`.
