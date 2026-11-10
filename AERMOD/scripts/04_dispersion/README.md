# AERMOD simulations

| Script | Source category |
|---|---|
| `run_drilling_completion.py` | Drilling and completion |
| `run_production.py` | Production leakage and heaters combined |
| `run_compressors.py` | Compressor facilities |
| `run_other_tier2.py` | Other Tier 2 facilities |

Each script reads the final inventories, the six parameter CSVs in `inputs/parameters`, AERMET meteorology, receptor DAT files and the AERMOD executable. `source_parameters.py` reads the source settings and writes the point/area source sections consistently across the four categories. Parameters are supplied in the tables rather than estimated during a model run.

```bash
python run.py run_compressors --config config/settings.json
```

Choose years, pollutants, source IDs and worker counts in the configuration. Use saved Python scripts for multiprocessing on HPC. Daily inventory fields are tonnes per day; monthly fields are total tonnes per month and are distributed evenly over the month's days. Point rates are g/s and area rates are g/s/m2; release weights sum to one. NO2 cases read the NOX inventory with the configured NO2 model options.

Results are stored by category/year/pollutant/source under `outputs/model_runs`. Validated receptor concentrations are compressed; detailed model output is removed after successful conversion and checks. Case status is consolidated in the category's records. Failed cases retain evidence for diagnosis. Successful existing cases are checked before being skipped.

The supplied case table defines the completed study cases. New sources or changed emissions require corresponding parameter/case entries. Before a full rerun with this portable code, compare representative point, area and multi-release inputs and results with the completed study cases. The new package has local syntax and synthetic checks, but has not been run through the complete HPC batch.
