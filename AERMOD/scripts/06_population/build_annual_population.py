"""Build annual total population on the fixed TriMAQs grid, using the supplied census archive.
Run this saved file through the accompanying SLURM job. Existing data are retained.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import sqlite3
import zipfile
import numpy as np
import pandas as pd
import geopandas as gpd
import population_spatial as spatial

# USER SETTINGS --------------------------------------------------------------
DATA_ROOT = Path(os.environ.get('TRIMAQS_DATA_ROOT', str(Path(__file__).resolve().parents[3]/'TriMAQs_data')))
PROJECT_ROOT = DATA_ROOT

CENSUS_ZIP = DATA_ROOT / 'inputs/source_data/census/Census_population.zip'
GRID_SPEC = DATA_ROOT / 'inputs/supporting/ml_grid/grid_spec_and_mask.npz'
GRID_CSV = DATA_ROOT / 'inputs/supporting/grids/study_grid_1km.csv'
OUTPUT_DIR = DATA_ROOT / 'inputs/population_demographics/general'
WORKERS = 20
# ---------------------------------------------------------------------------
VERSION = 'annual-total-population-v1'
require = spatial.require
sha256 = spatial.sha256
atomic_text = spatial.atomic_text


def extract_shape(archive, folder, stem, working):
    target = working / 'sources' / folder.replace('/', '_') / stem
    target.mkdir(parents=True, exist_ok=True)
    for ext in ['.shp', '.dbf', '.shx', '.prj', '.cpg']:
        matches = [n for n in archive.namelist() if n.replace('\\', '/').endswith('/' + folder + '/' + stem + ext)]
        if not matches and ext == '.cpg':
            continue
        require(len(matches) == 1, f'Missing or ambiguous source in ZIP: {folder}/{stem + ext}')
        info = archive.getinfo(matches[0])
        dest = target / (stem + ext)
        if dest.exists() and dest.stat().st_size == info.file_size:
            import zlib
            crc = 0
            with dest.open('rb') as f:
                for data in iter(lambda: f.read(8*1024*1024), b''):
                    crc = zlib.crc32(data, crc)
            if crc == info.CRC:
                continue
        temporary = dest.with_suffix(ext + '.partial')
        with archive.open(info) as source, temporary.open('wb') as output:
            for data in iter(lambda: source.read(8*1024*1024), b''):
                output.write(data)
        os.replace(temporary, dest)
    return target / (stem + '.shp')


def recover_recent_population(archive, frame, reference, working):
    """Recover total counts from 36 disjoint sex/age cells for ACS 2019-2023."""
    name = next(n for n in archive.namelist() if n.endswith('/Sex and age/nhgis0018_ts_nominal_tract.csv') and '/ACS ' in n)
    with archive.open(name) as stream:
        header = pd.read_csv(stream, nrows=0).columns
    columns = [c for c in header if re.fullmatch(r'B58[A-Z]{2}235', c)]
    require(len(columns) == 36, 'Unexpected sex/age schema; expected 36 mutually exclusive age/sex cells.')
    with archive.open(name) as stream:
        raw = pd.read_csv(stream, usecols=['GJOIN2023'] + columns, dtype={'GJOIN2023': str})
    raw['population'] = raw[columns].sum(axis=1, min_count=36)
    raw = raw.dropna(subset=['GJOIN2023'])
    require(not raw.GJOIN2023.duplicated().any(), 'Duplicate 2023 geography IDs in raw sex/age data.')
    lookup = raw.set_index('GJOIN2023').population
    # The source state selection omits some Maryland tracts. Use the supplied
    # 2020 reference polygons for missing study tracts, keeping this explicit.
    missing_ids = sorted(set(reference.GISJOIN) - set(frame.GISJOIN))
    added = reference[reference.GISJOIN.isin(missing_ids)][['GISJOIN', 'geometry']].copy()
    added['Total_Popu'] = added.GISJOIN.map(lookup)
    require(added.Total_Popu.notna().all(), 'Raw 2023 table cannot recover all missing reference tracts.')
    frame = pd.concat([frame, added], ignore_index=True)
    frame = gpd.GeoDataFrame(frame, geometry='geometry', crs=reference.crs)
    existing = pd.to_numeric(frame.Total_Popu, errors='coerce')
    recovered = frame.GISJOIN.map(lookup)
    check = frame.GISJOIN.isin(reference.GISJOIN) & existing.notna() & recovered.notna()
    require(np.allclose(existing[check], recovered[check], rtol=0, atol=0), '2023 total counts disagree with raw sex/age totals.')
    frame['Total_Popu'] = recovered.combine_first(existing)
    added[['GISJOIN', 'Total_Popu']].to_csv(working / 'recent_reference_geometry.csv', index=False)
    return frame, dict(added_reference_polygons=len(added), matching_existing_study_counts=int(check.sum())-len(added),
                       geometry_note='Missing study polygons use the supplied 2020 tract boundaries. This assumes those boundaries apply to these IDs in the 2019-2023 ACS period.')


def balanced_round(values, ids):
    require(np.isfinite(values).all() and (values >= 0).all(), 'Invalid annual population estimates.')
    result = np.floor(values).astype(np.int64)
    target = int(np.floor(float(values.sum()) + 0.5))
    remainder = target - int(result.sum())
    require(0 <= remainder <= len(result), 'Rounding remainder is inconsistent.')
    order = np.lexsort((ids, -(values - result)))
    result[order[:remainder]] += 1
    require(int(result.sum()) == target, 'Balanced rounding did not preserve the modeled total.')
    return result


def extrapolate(anchors, year):
    # Fit a line through the 2021 value to the 2017-2021 series. Restrict its
    # slope so predictions remain nonnegative throughout 2017-2023.
    last = anchors[2021]
    x = np.arange(-4, 1, dtype=float)
    recent = np.stack([anchors[y] for y in range(2017, 2022)])
    slope = (x[:, None] * (recent - last)).sum(axis=0) / np.square(x).sum()
    slope = np.clip(slope, -last / 2, last / 4)
    result = last + (year - 2021) * slope
    require((result >= -1e-10).all(), 'Constrained extrapolation became negative.')
    return np.maximum(result, 0)


def allocate(path, overlap, ids):
    frame = gpd.read_file(path)
    require('GISJOIN' in frame and 'Total_Popu' in frame, 'Missing population columns.')
    require(frame.GISJOIN.notna().all() and not frame.GISJOIN.duplicated().any(), 'Invalid source IDs.')
    values = pd.to_numeric(frame.Total_Popu, errors='coerce')
    lookup = pd.Series(values.to_numpy(), index=frame.GISJOIN.astype(str))
    links = pd.read_csv(overlap, dtype={'GISJOIN': str})
    population = links.GISJOIN.map(lookup)
    require(population.notna().all() and (population >= 0).all(), f'Missing/negative source counts within study grid: {path}')
    contributions = population.to_numpy() * links.area_weight.to_numpy()
    allocated = pd.Series(contributions).groupby(links.grid_ID).sum().reindex(ids, fill_value=0).to_numpy()
    return allocated


def build():
    require(CENSUS_ZIP.is_file(), f'Put Census_population.zip at {CENSUS_ZIP}, or change CENSUS_ZIP.')
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    working = OUTPUT_DIR / 'working'
    working.mkdir(exist_ok=True)
    # One build at a time; the operating system releases this lock on job exit.
    import fcntl
    lock = (working / 'build.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError('Another population build is using this output directory.')
    spatial.GRID_SPEC, spatial.GRID_CSV = GRID_SPEC, GRID_CSV
    spatial.OUTPUT_DIR = working / 'overlaps'
    spatial.OUTPUT_DIR.mkdir(exist_ok=True)
    spatial.WORKERS = WORKERS
    grid = spatial.load_grid()
    ids, xs, ys, lat, lon, grid_info = grid
    print(f'Population grid: {len(ids):,} cells; workers requested: {WORKERS}', flush=True)
    # ZIP identity uses size, member paths, CRCs and uncompressed sizes. Per-source
    # geometry hashes and raw-table CRCs are additionally retained below.
    with zipfile.ZipFile(CENSUS_ZIP) as archive:
        members = [(i.filename, i.CRC, i.file_size) for i in archive.infolist()]
        zip_identity = hashlib.sha256(json.dumps(members).encode()).hexdigest()
        identity = dict(version=VERSION, census_zip_identity=zip_identity, grid=grid_info)
        signature = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        identity_file = working / 'input_identity.json'
        if identity_file.exists():
            require(json.loads(identity_file.read_text())['signature'] == signature,
                    'Inputs changed. Choose a new OUTPUT_DIR to retain the previous build.')
        atomic_text(identity_file, json.dumps(dict(**identity, signature=signature), indent=2))
        db = sqlite3.connect(working / 'checkpoint.sqlite3')
        db.execute('CREATE TABLE IF NOT EXISTS chunks (signature TEXT,chunk INTEGER,payload BLOB,PRIMARY KEY(signature,chunk))')
        db.execute('CREATE TABLE IF NOT EXISTS finals (year INTEGER PRIMARY KEY,signature TEXT,sha256 TEXT)')
        db.execute('CREATE TABLE IF NOT EXISTS products (name TEXT PRIMARY KEY,sha256 TEXT)')
        ref_path = extract_shape(archive, 'Decennial', 'Total_pops_2020_tract', working)
        reference = gpd.read_file(ref_path).to_crs(5070)
        # Only add missing reference polygons that actually intersect study cells.
        from shapely.geometry import box
        from shapely.strtree import STRtree
        cells = [box(x-500,y-500,x+500,y+500) for x,y in zip(xs,ys)]
        tree = STRtree(cells)
        keep = []
        for geometry in reference.geometry:
            candidates = tree.query(geometry)
            keep.append(any(geometry.intersection(cells[int(i)] if isinstance(i,(int,np.integer)) else i).area > 1e-6 for i in candidates))
        reference = reference.loc[keep].copy()
        del cells, tree
        anchors, source_audit = {}, []
        for end in [2000] + list(range(2010, 2024)):
            model_year = 2000 if end == 2000 else end-2
            stem = 'Total_pops_2000_tract' if end == 2000 else f'Total_pops_{end-4}-{end}_tract'
            folder = 'Decennial' if end == 2000 else 'ACS 5-yr estimates'
            source = extract_shape(archive, folder, stem, working)
            note = {}
            if end == 2023:
                prepared = working / 'sources/recent_complete/Total_pops_2019-2023_tract.shp'
                if not (all(prepared.with_suffix(ext).is_file() for ext in ['.shp','.dbf','.shx','.prj'])
                        and (working / 'recent_source_note.json').is_file()):
                    prepared.parent.mkdir(exist_ok=True)
                    frame = gpd.read_file(source).to_crs(5070)[['GISJOIN','Total_Popu','geometry']]
                    frame, note = recover_recent_population(archive, frame, reference, working)
                    frame.to_file(prepared, index=False)
                    atomic_text(working / 'recent_source_note.json', json.dumps(note, indent=2))
                else:
                    note = json.loads((working / 'recent_source_note.json').read_text())
                source = prepared
            spatial.tract_file = lambda unused, source=source: source
            info = spatial.build_vintage(db, end, grid)
            require(info['grid_cells_without_tract_overlap'] == 0,
                    f'Source period ending {end} leaves study cells without tract coverage; stopped for review.')
            overlap = spatial.OUTPUT_DIR / info['output_file']
            cache = working / f'anchor_{model_year}.npy'
            old = db.execute('SELECT sha256 FROM products WHERE name=?',(cache.name,)).fetchone()
            if cache.exists():
                require(old is not None and old[0] == sha256(cache), 'Unverified anchor cache; choose a new OUTPUT_DIR.')
                values = np.load(cache, allow_pickle=False)
            else:
                values = allocate(source, overlap, ids)
                temporary = cache.with_suffix('.partial')
                with temporary.open('wb') as stream:
                    np.save(stream, values, allow_pickle=False)
                db.execute('INSERT OR REPLACE INTO products VALUES (?,?)',(cache.name,sha256(temporary)))
                db.commit()
                os.replace(temporary, cache)
            anchors[model_year] = values
            source_audit.append(dict(model_year=model_year, source_period='2000 census' if end==2000 else f'{end-4}-{end}',
                                     modeled_total=float(values.sum()), source_shp_sha256=sha256(source), **note))
            print(f'Anchor {model_year}: {values.sum():,.2f} allocated people', flush=True)
    general = OUTPUT_DIR / 'general'
    general.mkdir(exist_ok=True)
    summaries = []
    for year in range(2004,2024):
        if year < 2008:
            weight = (year-2000)/8
            values = (1-weight)*anchors[2000]+weight*anchors[2008]
            method = 'Linear interpolation between 2000 census and ACS 2006-2010 assigned to 2008'
        elif year <= 2021:
            values = anchors[year]
            method = 'ACS five-year estimate assigned to midpoint year'
        else:
            values = extrapolate(anchors, year)
            method = 'Constrained linear extrapolation through 2021 using 2017-2021 grid estimates'
        integer = balanced_round(values, ids)
        dest = general / f'population_{year}.csv'
        old = db.execute('SELECT sha256 FROM products WHERE name=?',(dest.name,)).fetchone()
        if dest.exists():
            require(old is not None and old[0] == sha256(dest), 'Existing population output cannot be verified; choose a new OUTPUT_DIR.')
        else:
            temporary = dest.with_suffix('.csv.partial')
            pd.DataFrame(dict(grid_ID=ids,latitude=lat,longitude=lon,population=integer,
                              population_estimate=values)).to_csv(temporary,index=False,float_format='%.15g')
            db.execute('INSERT OR REPLACE INTO products VALUES (?,?)',(dest.name,sha256(temporary)))
            db.commit()
            os.replace(temporary,dest)
        summaries.append(dict(year=year,method=method,grid_count=len(ids),modeled_total=float(values.sum()),
                              integer_total=int(integer.sum()),zero_cells=int((integer==0).sum()),
                              maximum_cell=int(integer.max()),file_sha256=sha256(dest)))
        print(f'{year}: {integer.sum():,} people; saved/verified {dest.name}',flush=True)
    pd.DataFrame(source_audit).to_csv(OUTPUT_DIR/'source_anchors.csv',index=False)
    pd.DataFrame(summaries).to_csv(OUTPUT_DIR/'annual_population_summary.csv',index=False)
    atomic_text(OUTPUT_DIR/'BUILD_SUMMARY.json',json.dumps(dict(state='complete_total_population',signature=signature,
        years=20,grid_cells=len(ids),nonnegative=True,integer_totals_preserved=True,
        limitations=['Annual totals are modeled, not calibrated to independent annual population controls.',
                     'ACS periods are assigned to midpoint years as a modeling convention.',
                     'Incomplete 2005-2009 ACS input is omitted.',
                     'Missing recent study polygons use supplied 2020 boundaries; IDs and counts are checked, but unchanged boundaries are assumed.',
                     'Demographics are not rebuilt by this total-population script.'],
        sources=source_audit,annual=summaries),indent=2))
    db.close()
    print(f'Complete: {general}. Share BUILD_SUMMARY.json and annual_population_summary.csv before replacing shared population inputs.',flush=True)
    lock.close()


if __name__ == '__main__':
    build()
