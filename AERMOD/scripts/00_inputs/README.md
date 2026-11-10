# Prepare the input package

`prepare_model_inputs.py` creates the final emission, parameter and supporting tables for the study data release. Readers who download those tables can proceed directly to meteorology and dispersion.

The stage reads completed model records and their AERMOD input files, the working inventory files referenced by those records, the compilation summary and the supporting-data archive. Set `export.source_model_runs`, `export.source_archive`, `export.compilation_summary`, `export.destination_inputs` and `export.report_folder` in the configuration.

```bash
python run.py prepare_model_inputs --config config/settings.json
```

The destination contains `emissions`, `parameters` and `supporting`, with folder READMEs, a data dictionary and a file manifest. Daily emissions are tonnes per day; every monthly inventory is converted to tonnes per calendar month. Annual totals are checked against the completed model records. Source locations, release dimensions and grid selections are taken from the final modeled settings. Sources outside the study boundary are excluded.

The six parameter CSVs are described in `docs/INPUTS.md`. These are ordinary model inputs, read directly by the four dispersion scripts. No supplementary JSON parameter datasets are needed.

The script requires a new destination and leaves existing data unchanged. It checks source counts against the compilation summary and stops on conflicting records, changed inventories or unmatched mass totals. A single build summary is written outside the shared inputs. Review that summary before publishing. Existing uncertainties in source assignments and estimated release parameters remain relevant and are recorded in the table notes.
