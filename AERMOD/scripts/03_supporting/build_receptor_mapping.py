"""Build the four-receptor sampling map for each study grid."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
from pathlib import Path
from datetime import datetime
import json, hashlib
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from pyproj import Transformer
BASE = Path(_configured_path('build_receptor_mapping.base_1'))
GRID_MAPPING = Path(_configured_path('grid_links'))
GRID = BASE / 'inputs/supporting/grids/study_grid_1km.csv'
AERMAP = BASE / 'inputs/supporting/receptors'
OUTPUT = BASE / 'inputs/supporting/receptor_mapping'
XY = ['utm_easting', 'utm_northing']

def require(ok, message):
    if not ok:
        raise ValueError(message)

def file_record(path):
    info = path.stat()
    return dict(path=str(path), size=info.st_size, mtime_ns=info.st_mtime_ns, sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def read_receptors(path):
    points = []
    for line in path.read_text().splitlines():
        tokens = line.split()
        if tokens and tokens[0] == 'RE':
            tokens = tokens[1:]
        if tokens and tokens[0] in {'EVALCART', 'DISCCART'}:
            points.append((round(float(tokens[1]), 2), round(float(tokens[2]), 2)))
    values = np.unique(np.asarray(points, dtype=float).reshape(-1, 2), axis=0)
    require(len(values) >= 4 and np.isfinite(values).all(), 'Invalid/fewer than four receptors: ' + str(path))
    return values

def main():
    require(not OUTPUT.exists(), 'Output already exists; retain it and choose a new OUTPUT if rebuilding.')
    None
    old = pd.read_csv(GRID_MAPPING)
    links = old[['grid_ID_1km', 'grid_ID_10km']].drop_duplicates()
    require(not links.grid_ID_1km.duplicated().any(), 'Conflicting 1km-to-10km links')
    grid = pd.read_csv(GRID).rename(columns={'grid_ID': 'grid_ID_1km'})
    require(not grid.grid_ID_1km.duplicated().any(), 'Duplicate 1km grid IDs')
    grid = grid[['grid_ID_1km', 'lat', 'lon']].merge(links, on='grid_ID_1km', how='left', validate='one_to_one')
    require(grid.grid_ID_10km.notna().all(), 'Some output grids lack a parent 10km grid')
    transform = Transformer.from_crs('EPSG:4326', 'EPSG:32618', always_xy=True)
    grid['center_easting'], grid['center_northing'] = transform.transform(grid.lon.to_numpy(), grid.lat.to_numpy())
    frames = []
    reference = []
    inputs = []
    for position, (gid, group) in enumerate(grid.groupby('grid_ID_10km'), 1):
        gid = int(gid)
        path = AERMAP / f'Grid_{gid}' / f'receptor_output_{gid}_grid.DAT'
        require(path.is_file(), 'Missing current AERMAP file: ' + str(path))
        points = read_receptors(path)
        reference.append(pd.DataFrame(points, columns=XY))
        inputs.append(file_record(path))
        distance, indices = cKDTree(points).query(group[['center_easting', 'center_northing']].to_numpy(), k=4)
        frame = group.loc[group.index.repeat(4), ['grid_ID_1km', 'grid_ID_10km', 'lat', 'lon']].reset_index(drop=True)
        coords = points[indices.ravel()]
        frame['utm_easting'] = coords[:, 0]
        frame['utm_northing'] = coords[:, 1]
        frame['distance_to_grid_center_m'] = distance.ravel()
        frames.append(frame)
        if position % 100 == 0:
            print('AERMAP grids read:', position, flush=True)
    selected = pd.concat(frames, ignore_index=True)
    reference = pd.concat(reference, ignore_index=True).drop_duplicates().reset_index(drop=True)
    require((selected.groupby('grid_ID_1km').size() == 4).all(), 'Not four assignments per output grid')
    require(not selected.duplicated(XY + ['grid_ID_1km']).any(), 'Duplicate receptor within a grid')
    index = pd.MultiIndex.from_frame(reference)
    checks = []
    # Every model CSV is checked against this map during compilation.
    report = dict(state='ready' if all((c['unmatched_coordinates'] == 0 for c in checks)) else 'blocked', created=datetime.now().isoformat(), grid_count=len(grid), assignment_count=len(selected), unique_reference_receptors=len(reference), source_checks=checks, sampling_distance_m=selected.distance_to_grid_center_m.describe(percentiles=[0.5, 0.95, 0.99]).to_dict(), inputs=[file_record(GRID_MAPPING), file_record(GRID)] + inputs, assumptions='Four nearest current receptors within the inherited parent 10km receptor file, equal weights. Receptor reuse among target cells is intentional. This estimates center-neighborhood concentrations, not exact cell area means. Ties follow cKDTree ordering.', limits='Checks all coordinates of supplied failing-case examples, not every model CSV. Compiler still checks every CSV during compilation. AERMOD outputs and old mapping files remain unchanged.')
    OUTPUT.mkdir(parents=True, exist_ok=False)
    if report['state'] == 'ready':
        selected.to_csv(OUTPUT / 'receptor_mapping.csv', index=False)
        reference.to_csv(OUTPUT / 'reference_receptors.csv', index=False)
        selected.groupby('grid_ID_1km').agg(receptors=('utm_easting', 'size'), maximum_distance_m=('distance_to_grid_center_m', 'max')).to_csv(OUTPUT / 'grid_sampling_summary.csv')
    report_path=BASE/'working/receptor_mapping_summary.json'
    report_path.parent.mkdir(parents=True,exist_ok=True)
    report_path.write_text(json.dumps(report,indent=2))
    print('Report:',report_path,flush=True)
    require(report['state'] == 'ready', 'Unmatched current coordinates remain; share report. No mapping CSVs published.')
    print('Mappings saved. Share BUILD_REPORT.json before recompilation; existing compiled groups used the older mapping.', flush=True)
if __name__ == '__main__':
    main()
