# Terrain, receptors and spatial mapping

`run_aermap.py` calculates terrain information at the supplied receptor coordinates using elevation data and the AERMAP executable. It reads the coordinate templates from `run_aermap.coordinate_templates`; it does not define a different receptor grid. The output must be separate from those templates. Elevation rasters must share EPSG:4269 for this stage. Existing receptor DAT files can be used directly when running AERMOD, without repeating the terrain calculation.

`build_receptor_mapping.py` selects four current receptors for each 1 km study cell, using the cell's associated 10 km receptor grid. The study grid CSV supplies the cell centers. The 1km-to-10km link table is set by `grid_links` in the configuration.

The mapping outputs are `receptor_mapping.csv` and `reference_receptors.csv` under `inputs/supporting/receptor_mapping`; processing summaries are kept separately. One receptor may serve several study cells. Concentrations must be summed at unique receptors before averaging the selected receptors for a grid.

```bash
python run.py run_aermap --config config/settings.json
python run.py build_receptor_mapping --config config/settings.json
```

The selected receptor mean is a spatial estimate, not an exact area-weighted cell mean. The final NetCDF raster uses the collaborator grid in `inputs/supporting/ml_grid/grid_spec_and_mask.npz`: x/y are indices, lat/lon are geographic coordinates, and mask=True is outside the study area. The 1 km grid_ID is the zero-based flattened raster index.
