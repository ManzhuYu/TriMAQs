"""Recreate tract-to-study-grid overlap tables. Run this saved file on HPC.

This step calculates areas only; it does not interpolate or round population.
Existing census files and population estimates are never changed.
"""
from pathlib import Path
import concurrent.futures
import csv
import hashlib
import io
import json
import os
import sqlite3
import time
import zipfile
import zlib

import numpy as np
import pandas as pd
import geopandas as gpd
from pyproj import Transformer
from shapely import wkb
from shapely.geometry import box
from shapely.strtree import STRtree
from shapely.validation import make_valid

# USER SETTINGS -------------------------------------------------------------
DATA_ROOT = Path(os.environ.get('TRIMAQS_DATA_ROOT', str(Path(__file__).resolve().parents[3]/'TriMAQs_data')))
PROJECT_ROOT = DATA_ROOT

GRID_SPEC = DATA_ROOT / 'inputs/supporting/ml_grid/grid_spec_and_mask.npz'
GRID_CSV = DATA_ROOT / 'inputs/supporting/grids/study_grid_1km.csv'
CENSUS_ROOT = PROJECT_ROOT / 'input_data/Population/census tract'
# Used only if the decennial shapefiles are not already in CENSUS_ROOT.
CENSUS_ZIP = DATA_ROOT / 'inputs/source_data/census/Census_population.zip'
OUTPUT_DIR = PROJECT_ROOT / 'input_data/Population/Population_overlaps'
YEARS = [2000, 2010, 2020]
WORKERS = 20
CHUNK_SIZE = 2500
LATTICE_TOLERANCE_M = 2.5
# ---------------------------------------------------------------------------
VERSION = 'tract-overlap-1km-v1'
COLUMNS = ['grid_ID', 'GISJOIN', 'overlap_area_km2', 'tract_area_km2', 'area_weight']
TRACTS = TREE = GEOM_INDEX = None


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def atomic_text(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.partial')
    temporary.write_text(value, encoding='utf-8')
    os.replace(temporary, path)


def tract_file(year):
    name = f'Total_pops_{year}_tract'
    candidates = [CENSUS_ROOT / 'Decennial' / (name + '.shp'),
                  CENSUS_ROOT / 'census tract/Decennial' / (name + '.shp')]
    for candidate in candidates:
        if all(candidate.with_suffix(ext).is_file() for ext in ['.shp', '.dbf', '.shx', '.prj']):
            return candidate
    require(CENSUS_ZIP.is_file(), f'Cannot find {name}.shp. Set CENSUS_ROOT or CENSUS_ZIP at the top of this file.')
    destination = OUTPUT_DIR / 'working/tract_sources' / str(year)
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(CENSUS_ZIP) as archive:
        for ext in ['.shp', '.dbf', '.shx', '.prj', '.cpg']:
            matches = [n for n in archive.namelist()
                       if n.replace('\\', '/').endswith('/Decennial/' + name + ext)]
            if not matches and ext == '.cpg':
                continue
            require(len(matches) == 1, f'Expected exactly one {name + ext} in the ZIP.')
            target = destination / (name + ext)
            info = archive.getinfo(matches[0])
            if target.is_file() and target.stat().st_size == info.file_size:
                with target.open('rb') as stream:
                    crc = 0
                    for data in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                        crc = zlib.crc32(data, crc)
                if crc == info.CRC:
                    continue
            temporary = target.with_suffix(ext + '.partial')
            with archive.open(info) as source, temporary.open('wb') as stream:
                for data in iter(lambda: source.read(8 * 1024 * 1024), b''):
                    stream.write(data)
            os.replace(temporary, target)
    return destination / (name + '.shp')


def load_grid():
    require(GRID_SPEC.is_file(), f'Missing grid specification: {GRID_SPEC}')
    with np.load(GRID_SPEC, allow_pickle=False) as data:
        lat, lon = data['lat'], data['lon']
        mask = data['study_mask'].astype(bool)
    require(lat.ndim == 2 and lat.shape == lon.shape == mask.shape, 'Grid arrays must have matching two-dimensional shapes.')
    require(np.isfinite(lat).all() and np.isfinite(lon).all(), 'Grid coordinates contain missing values.')
    x, y = Transformer.from_crs(4326, 5070, always_xy=True).transform(lon, lat)
    rows, cols = np.indices(lat.shape)
    sx, sy = int(np.sign(np.median(np.diff(x, axis=1)))), int(np.sign(np.median(np.diff(y, axis=0))))
    require(sx != 0 and sy != 0, 'Cannot determine grid orientation.')
    ox, oy = float(np.median(x - cols * sx * 1000)), float(np.median(y - rows * sy * 1000))
    fx, fy = ox + cols * sx * 1000, oy + rows * sy * 1000
    error = np.hypot(fx - x, fy - y)
    require(float(error.max()) <= LATTICE_TOLERANCE_M,
            f'Grid is not the expected regular 1 km lattice: maximum residual {error.max():.3f} m.')
    inside = ~mask
    ids = (rows * lat.shape[1] + cols)[inside].astype('int64')
    require(len(ids) > 0, 'Study mask selects no cells.')
    if GRID_CSV.is_file():
        frame = pd.read_csv(GRID_CSV)
        lookup = {str(c).lower(): c for c in frame.columns}
        require('grid_id' in lookup, 'Grid CSV has no grid_ID column.')
        csv_ids = pd.to_numeric(frame[lookup['grid_id']], errors='raise').to_numpy()
        require(np.isfinite(csv_ids).all() and np.equal(csv_ids, np.floor(csv_ids)).all(), 'Noninteger grid IDs.')
        require(len(csv_ids) == len(ids) and len(np.unique(csv_ids)) == len(ids)
                and np.array_equal(np.sort(csv_ids), ids), 'Grid CSV IDs do not match the NPZ study mask.')
        ordered = frame.assign(_id=csv_ids).set_index('_id').loc[ids]
        for keys, expected in [(('lat', 'latitude'), lat[inside]), (('lon', 'longitude'), lon[inside])]:
            column = next((lookup[k] for k in keys if k in lookup), None)
            if column:
                require(np.allclose(pd.to_numeric(ordered[column]), expected, atol=5e-6, rtol=0),
                        f'Grid CSV {column} does not match the NPZ.')
    audit = dict(crs='EPSG:5070', cell_width_m=1000, grid_count=len(ids),
                 origin_x_m=ox, origin_y_m=oy, x_direction=sx, y_direction=sy,
                 max_center_residual_m=float(error.max()), rms_center_residual_m=float(np.sqrt(np.mean(error**2))),
                 grid_spec_sha256=sha256(GRID_SPEC),
                 grid_csv_sha256=sha256(GRID_CSV) if GRID_CSV.is_file() else None,
                 mask_convention='True outside study area; grid_ID = row * number_of_columns + column')
    return ids, fx[inside], fy[inside], lat[inside], lon[inside], audit


def load_tracts(path):
    frame = gpd.read_file(path)
    require(frame.crs is not None, f'Missing tract CRS: {path}')
    require('GISJOIN' in frame, f'Missing GISJOIN: {path}')
    require(frame.GISJOIN.notna().all(), 'Tract GISJOIN contains missing values.')
    frame['GISJOIN'] = frame.GISJOIN.astype(str).str.strip()
    require(frame.GISJOIN.ne('').all() and not frame.GISJOIN.duplicated().any(), 'Empty or duplicate tract GISJOIN.')
    frame = frame.to_crs(5070).sort_values('GISJOIN').reset_index(drop=True)
    require(frame.geometry.notna().all() and (~frame.geometry.is_empty).all(), 'Missing tract geometry.')
    require(frame.geom_type.isin(['Polygon', 'MultiPolygon']).all(), 'Unexpected tract geometry type.')
    corrections = []
    for index in frame.index[~frame.geometry.is_valid]:
        original = frame.at[index, 'geometry']
        geometry = make_valid(original)
        delta = abs(geometry.area - original.area)
        require(geometry.geom_type in ['Polygon', 'MultiPolygon'] and geometry.is_valid
                and delta <= max(1.0, original.area * 1e-8),
                f'Tract {frame.at[index, "GISJOIN"]} requires a substantial geometry change; stopped for review.')
        corrections.append(dict(GISJOIN=frame.at[index, 'GISJOIN'], area_change_m2=delta))
        frame.at[index, 'geometry'] = geometry
    atomic_text(OUTPUT_DIR / ('geometry_processing_' + path.stem + '.json'), json.dumps(corrections, indent=2))
    areas = frame.geometry.area.to_numpy() / 1e6
    require(np.isfinite(areas).all() and (areas > 0).all(), 'Invalid tract areas.')
    return frame.GISJOIN.tolist(), [g.wkb for g in frame.geometry], areas


def initialize_worker(geometries):
    global TRACTS, TREE, GEOM_INDEX
    TRACTS = [wkb.loads(g) for g in geometries]
    TREE = STRtree(TRACTS)
    GEOM_INDEX = {id(g): i for i, g in enumerate(TRACTS)}


def calculate_chunk(task):
    chunk, ids, xs, ys = task
    records = []
    for grid_id, x, y in zip(ids, xs, ys):
        cell = box(x - 500, y - 500, x + 500, y + 500)
        candidates = TREE.query(cell)
        indexes = [int(item) if isinstance(item, (int, np.integer)) else GEOM_INDEX[id(item)] for item in candidates]
        for index in sorted(indexes):
            area = cell.intersection(TRACTS[index]).area / 1e6
            if area > 0:
                records.append((int(grid_id), index, area))
    return chunk, zlib.compress(json.dumps(records, separators=(',', ':')).encode(), level=3)


def build_vintage(database, year, grid):
    ids, xs, ys, lat, lon, grid_audit = grid
    path = tract_file(year)
    tract_ids, geometries, areas = load_tracts(path)
    sources = {ext: sha256(path.with_suffix(ext)) for ext in ['.shp', '.dbf', '.shx', '.prj']}
    provenance = dict(version=VERSION, year=year, grid=grid_audit, sources=sources, chunk_size=CHUNK_SIZE)
    signature = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    filename = 'Overlap_Areas_master.csv' if year == 2020 else f'Overlap_Areas_{year}.csv'
    destination = OUTPUT_DIR / filename
    old = database.execute('SELECT signature, sha256 FROM finals WHERE year=?', (year,)).fetchone()
    if destination.exists():
        require(old is not None and old[0] == signature and old[1] == sha256(destination),
                f'Existing output cannot be verified: {destination}. Choose a new OUTPUT_DIR; originals retained.')
        print(f'{year}: verified existing {filename}', flush=True)
        return json.loads((OUTPUT_DIR / f'summary_{year}.json').read_text())
    cached = {row[0] for row in database.execute('SELECT chunk FROM chunks WHERE signature=?', (signature,))}
    tasks = [(i // CHUNK_SIZE, ids[i:i+CHUNK_SIZE], xs[i:i+CHUNK_SIZE], ys[i:i+CHUNK_SIZE])
             for i in range(0, len(ids), CHUNK_SIZE) if i // CHUNK_SIZE not in cached]
    total_chunks = (len(ids) + CHUNK_SIZE - 1) // CHUNK_SIZE
    workers = min(WORKERS, int(os.environ.get('SLURM_CPUS_PER_TASK', WORKERS)))
    require(workers >= 1, 'WORKERS must be positive.')
    print(f'{year}: {len(tract_ids):,} tracts; {len(cached)}/{total_chunks} chunks reused; workers={workers}', flush=True)
    def save(result):
        chunk, payload = result
        database.execute('INSERT INTO chunks VALUES (?,?,?)', (signature, chunk, payload))
        database.commit()
        print(f'{year}: saved chunk {chunk+1}/{total_chunks}', flush=True)
    if tasks and workers == 1:
        initialize_worker(geometries)
        for task in tasks:
            save(calculate_chunk(task))
    elif tasks:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers, initializer=initialize_worker,
                                                    initargs=(geometries,)) as pool:
            futures = [pool.submit(calculate_chunk, task) for task in tasks]
            for future in concurrent.futures.as_completed(futures):
                save(future.result())
    covered_grid = np.zeros(len(ids), dtype=float)
    covered_tract = np.zeros(len(tract_ids), dtype=float)
    temporary = destination.with_name(destination.name + '.partial')
    row_count = 0
    with temporary.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(COLUMNS)
        for chunk in range(total_chunks):
            value = database.execute('SELECT payload FROM chunks WHERE signature=? AND chunk=?', (signature, chunk)).fetchone()
            require(value is not None, f'Missing checkpoint chunk {chunk}.')
            for grid_id, index, area in json.loads(zlib.decompress(value[0])):
                position = int(np.searchsorted(ids, grid_id))
                require(position < len(ids) and ids[position] == grid_id, 'Checkpoint grid ID mismatch.')
                covered_grid[position] += area
                covered_tract[index] += area
                writer.writerow([grid_id, tract_ids[index], format(area, '.12g'), format(areas[index], '.12g'), format(area/areas[index], '.12g')])
                row_count += 1
    # 1 m2 allowance for numerical geometry operations; never renormalize areas.
    require((covered_grid <= 1 + 1e-6).all(), 'Overlapping tracts overfill one or more grid cells; final CSV not published.')
    require((covered_tract <= areas + 1e-6).all(), 'Overlapping grid cells overfill a tract; final CSV not published.')
    pd.DataFrame(dict(grid_ID=ids, latitude=lat, longitude=lon, covered_area_km2=covered_grid,
                      grid_area_km2=np.ones(len(ids)))).to_csv(OUTPUT_DIR / f'grid_coverage_{year}.csv', index=False)
    pd.DataFrame(dict(GISJOIN=tract_ids, tract_area_km2=areas, covered_area_km2=covered_tract,
                      coverage_fraction=covered_tract/areas)).to_csv(OUTPUT_DIR / f'tract_coverage_{year}.csv', index=False)
    summary = dict(**provenance, signature=signature, output_file=filename, rows=row_count,
                   grid_cells_without_tract_overlap=int((covered_grid == 0).sum()),
                   partially_covered_grid_cells=int(((covered_grid > 0) & (covered_grid < 1-1e-6)).sum()),
                   covered_area_km2=float(covered_grid.sum()), output_sha256=sha256(temporary),
                   state='complete_geometry_only', note='Use matching tract vintage. Full tract area is the population-allocation denominator; do not renormalize partial coverage.')
    atomic_text(OUTPUT_DIR / f'summary_{year}.json', json.dumps(summary, indent=2))
    # Commit expected final identity before rename, so a stopped job can recover.
    database.execute('INSERT OR REPLACE INTO finals VALUES (?,?,?)', (year, signature, summary['output_sha256']))
    database.commit()
    os.replace(temporary, destination)
    print(f'{year}: saved {filename}; {row_count:,} rows; {summary["grid_cells_without_tract_overlap"]} cells without tract coverage', flush=True)
    return summary


def main():
    require(CHUNK_SIZE > 0 and WORKERS > 0, 'Invalid worker or chunk settings.')
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    working = OUTPUT_DIR / 'working'
    working.mkdir(exist_ok=True)
    # An advisory lock is released automatically even when a SLURM job stops.
    # Run on Linux/HPC; Windows tests call the functions directly.
    import fcntl
    with (working / 'build.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another overlap build is using this OUTPUT_DIR.')
        grid = load_grid()
        atomic_text(OUTPUT_DIR / 'grid_geometry.json', json.dumps(grid[-1], indent=2))
        database = sqlite3.connect(working / 'overlap_checkpoint.sqlite3')
        database.execute('CREATE TABLE IF NOT EXISTS chunks (signature TEXT, chunk INTEGER, payload BLOB, PRIMARY KEY(signature,chunk))')
        database.execute('CREATE TABLE IF NOT EXISTS finals (year INTEGER PRIMARY KEY, signature TEXT, sha256 TEXT)')
        summaries = [build_vintage(database, year, grid) for year in YEARS]
        database.close()
        atomic_text(OUTPUT_DIR / 'BUILD_SUMMARY.json', json.dumps(summaries, indent=2))
        print(f'Complete. Share BUILD_SUMMARY.json and the three grid_coverage CSV files from {OUTPUT_DIR}', flush=True)
        print('Population estimates have not been rebuilt. No source files were changed.', flush=True)


if __name__ == '__main__':
    main()
