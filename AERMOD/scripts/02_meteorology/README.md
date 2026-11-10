# Meteorological inputs

These scripts prepare ERA5 data and run AERMET for each study year and meteorological grid. They are needed when rebuilding meteorology; prepared SFC/PFL files can be used directly by the dispersion workflow.

| Script | Purpose |
|---|---|
| `prepare_era5.py` | Read raw surface and pressure-level products and write annual inputs |
| `inspect_surface_fields.py` | Identify missing or invalid surface values |
| `complete_surface_fields.py` | Apply the documented small-gap treatment and record substitutions |
| `run_aermet.py` | Run Stage 1 and Stage 2 for the selected grids and years |
| `check_aermet.py` | Check that the expected year/grid results are present and usable |

Input locations, executables, years and worker counts are in the shared configuration. Raw products and prepared products are stored separately. ERA5 timestamps remain UTC during preprocessing; the model input convention uses fixed local standard time. Neighboring-year records are needed at annual boundaries.

The surface completion step uses its inspection report. Set the report path to the file created by `inspect_surface_fields`. The study rules include zeroing very small negative precipitation increments within the stated tolerance, same-hour spatial donors for identified radiation gaps within 30 km, and nearest-time values for the two specified missing boundary hours. These estimates are recorded; they should not be described as retrieved ERA5 observations.

AERMET uses pressure-level profiles, terrain elevation and AERSURFACE inputs. Below-terrain levels and pressure values that become duplicates after formatting are handled explicitly. The output check distinguishes documented missing-hour cases from complete meteorological records. Do not fill missing meteorology by assuming zero concentrations.

```bash
python run.py prepare_era5 --config config/settings.json
python run.py inspect_surface_fields --config config/settings.json
python run.py complete_surface_fields --config config/settings.json
python run.py run_aermet --config config/settings.json
python run.py check_aermet --config config/settings.json
```

Full meteorological files are large and may be supplied separately from the main dataset. Model runs also require their associated AERMET completion and quality records, not only the SFC/PFL filenames.
