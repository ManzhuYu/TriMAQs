"""Compile daily and annual grid concentrations for each source category."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
from pathlib import Path
from datetime import datetime
from decimal import Decimal
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
import os, json, csv, sqlite3, calendar, re, tempfile, shutil
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from pyproj import Transformer
BASE = Path(_configured_path('compile_concentrations.base_1'))
INPUT = BASE / 'inputs'
ROOT = BASE / 'outputs/model_runs'
MAPPING = BASE / 'inputs/supporting/receptor_mapping/receptor_mapping.csv'
REFERENCE_MAPPING = BASE / 'inputs/supporting/receptor_mapping/reference_receptors.csv'
GRID = INPUT / 'supporting/grids/study_grid_1km.csv'
BOUNDARY = INPUT / 'supporting/boundary/study_boundary.shp'
OUTPUT = BASE / 'outputs/compiled'
SCRATCH = Path(_configured_path('compile_concentrations.scratch_2'))
MAX_WORKERS = _setting('compile_concentrations.max_workers')
RETRY_SUMMARY = None
YEARS = _setting('compile_concentrations.years')
POLLUTANTS = _setting('compile_concentrations.pollutants')
WORKFLOWS = ['drilling_completion', 'production', 'compressors', 'other_tier2']
ENV = None

def require(ok, message):
    if not ok:
        raise ValueError(message)

def ident(value):
    text = str(value).strip()
    try:
        return str(int(Decimal(text)))
    except Exception:
        return re.sub('\\D', '', text) if re.fullmatch('[\\d-]+', text) else text

def inventory_path(workflow,year,pollutant):
    p='NOX' if pollutant in {'NOX','NO2'} else pollutant
    token='pm25' if p=='PM2.5' else p.lower()
    prefix='compressor' if workflow=='compressors' else workflow
    return INPUT/f'emissions/{workflow}/{p}/{prefix}_emissions_{token}_{year}.csv'

def load_environment():
    mapping = pd.read_csv(MAPPING, dtype=str)
    if 'grid_ID_1km' not in mapping and 'grid_ID' in mapping:
        mapping = mapping.rename(columns={'grid_ID': 'grid_ID_1km'})
    require({'utm_easting', 'utm_northing', 'grid_ID_1km'}.issubset(mapping), 'Mapping requires utm_easting, utm_northing, and grid_ID_1km or grid_ID')
    require(mapping['grid_ID_1km'].notna().all(), 'Missing receptor grid assignments')
    for col in ['utm_easting', 'utm_northing']:
        mapping[col] = pd.to_numeric(mapping[col], errors='raise').round(2)
    require(np.isfinite(mapping[['utm_easting', 'utm_northing']].to_numpy()).all(), 'Invalid mapping coordinates')
    mapping['grid_ID_1km'] = mapping.grid_ID_1km.map(ident)
    mapping = mapping[['utm_easting', 'utm_northing', 'grid_ID_1km']].drop_duplicates()
    reference = pd.read_csv(REFERENCE_MAPPING, usecols=['utm_easting', 'utm_northing'])
    for col in ['utm_easting', 'utm_northing']:
        reference[col] = pd.to_numeric(reference[col], errors='raise').round(2)
    require(np.isfinite(reference.to_numpy()).all(), 'Invalid reference receptor coordinates')
    reference = pd.concat([reference, mapping[['utm_easting', 'utm_northing']]], ignore_index=True).drop_duplicates()
    index = pd.MultiIndex.from_frame(reference)
    positions = index.get_indexer(pd.MultiIndex.from_frame(mapping[['utm_easting', 'utm_northing']]))
    grids = pd.read_csv(GRID, dtype=str)
    grids.columns = grids.columns.str.lower().str.strip()
    key = 'grid_id' if 'grid_id' in grids else 'grid_id_1km'
    gids = grids[key].map(ident).tolist()
    require(len(set(gids)) == len(gids), 'Duplicate output grid IDs')
    lookup = {g: i for i, g in enumerate(gids)}
    codes = np.asarray([lookup.get(g, -1) for g in mapping.grid_ID_1km], dtype=np.int64)
    counts = np.bincount(codes[codes >= 0], minlength=len(gids))
    county = gpd.read_file(BOUNDARY)
    require(county.crs is not None and county.geometry.is_valid.all(), 'Invalid study boundary/CRS')
    county = county.to_crs('EPSG:4326')
    domain = county.geometry.union_all() if hasattr(county.geometry, 'union_all') else county.geometry.unary_union
    source_review={}
    names={'drilling_completion':'drilling_completion_data.csv','production':'production_data.csv','compressors':'compressor_data_with_elevation.csv','other_tier2':'other_tier2_data_with_elevation.csv'}
    for category,name in names.items():
        frame=pd.read_csv(INPUT/'parameters'/name,dtype=str).fillna('')
        source_review[category]=set(frame.loc[frame.location_note.ne(''),'source_id'])
    tables = {}
    return dict(source_review=source_review,index=index, positions=positions, codes=codes, counts=counts, gids=gids, domain=domain, tables=tables, transformer=Transformer.from_crs('EPSG:32618', 'EPSG:4326', always_xy=True))

def initialize():
    global ENV
    ENV = load_environment()

def source_points(workflow,source,row,data):
    if data and data.get('source_release_parameters'):
        return [(float(r['longitude']),float(r['latitude'])) for r in data['source_release_parameters']]
    lat=row.get('latitude_decimal',row.get('latitude',row.get('lat')))
    lon=row.get('longitude_decimal',row.get('longitude',row.get('lon')))
    require(lat is not None and lon is not None,'Source position missing: '+source)
    return [(float(lon),float(lat))]

def read_done(folder, db):
    loose = folder / 'AERMOD_DONE.json'
    if loose.is_file():
        return json.loads(loose.read_text())
    entry = db.execute('SELECT record FROM cases WHERE source=?', (folder.name,)).fetchone() if db else None
    return json.loads(entry[0]) if entry else None

def emission_columns(frame,workflow,year):
    if workflow=='drilling_completion':return [c for c in frame if re.fullmatch(str(year)+r'-\d{2}-\d{2}',c)]
    return [calendar.month_abbr[m].lower()+'_tonnes' for m in range(1,13) if calendar.month_abbr[m].lower()+'_tonnes' in frame]

def compile_group(task):
    workflow, year, pollutant = task
    config = dict(workflow=workflow, year=year, pollutant=pollutant)
    source_inventory = inventory_path(*task)
    if not source_inventory.is_file():
        return dict(**config, state='inventory_absent', path=str(source_inventory))
    days = 366 if calendar.isleap(year) else 365
    columns = [f'day{i}_conc' for i in range(1, days + 1)]
    destination = OUTPUT / workflow / pollutant
    stem = f'{year}_{pollutant}'
    daily = destination / (stem + '_daily.csv.gz')
    annual = destination / (stem + '_annual.csv')
    audit = destination / (stem + '_audit.json')
    require(not any((p.exists() for p in [daily, annual, audit])), 'Compiled outputs exist; select new OUTPUT or remove only an explicitly rejected compilation')
    frame = pd.read_csv(source_inventory, dtype=str)
    frame.columns = frame.columns.str.lower().str.strip()
    key = 'source_id'
    ids = frame[key].map(ident)
    if workflow == 'drilling_completion':
        ids = ids.str.zfill(10)
    require(not ids.duplicated().any(), 'Duplicate source IDs in ' + str(source_inventory))
    emissions = emission_columns(frame, workflow, year)
    require(emissions, 'No recognized emission columns')
    numbers = frame[emissions].apply(pd.to_numeric, errors='raise').fillna(0.0)
    require(np.isfinite(numbers.to_numpy()).all() and numbers.ge(0).all().all(), 'Invalid emission entries')
    positive = numbers.sum(axis=1).gt(0)
    group = ROOT / workflow / str(year) / pollutant
    db = None
    if (group / 'audit.sqlite3').is_file():
        db = sqlite3.connect((group / 'audit.sqlite3').as_uri() + '?mode=ro', uri=True, timeout=120)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f'{workflow}_{year}_{pollutant}_', dir=SCRATCH))
    accumulated = None
    included = 0
    excluded = []
    met_review = 0
    points_processed = 0
    source_review = []
    try:
        accumulated = np.memmap(temporary / 'receptors.f64', dtype='float64', mode='w+', shape=(len(ENV['index']), days))
        accumulated[:] = 0
        for position in np.flatnonzero(positive.to_numpy()):
            source = ids.iloc[position]
            row = frame.iloc[position].to_dict()
            folder = group / source
            data = read_done(folder, db)
            require(data is not None and data.get('state') == 'completed', 'Missing successful case: ' + str(folder))
            locations = source_points(workflow, source, row, data)
            require(all((np.isfinite([lon, lat]).all() and -180 <= lon <= 180 and (-90 <= lat <= 90) for lon, lat in locations)), 'Invalid source coordinates')
            inside = [ENV['domain'].covers(Point(lon, lat)) for lon, lat in locations]
            if not any(inside):
                excluded.append(source)
                continue
            require(all(inside), 'Mixed inside/outside release origins need explicit emission allocation review: ' + source)
            require(data is not None and data.get('state') == 'completed', 'Missing successful case: ' + str(folder))
            require(data.get('folder') == str(folder), 'Completion record path differs: ' + str(folder))
            require(str(data['source_id']) == source and int(data['year']) == year and (data['pollutant'] == pollutant), 'Case identity mismatch')
            require(data['model']['fatal_errors'] == 0, 'AERMOD fatal errors recorded')
            unit_label = str(data['model'].get('units', '')).upper()
            require('MICROGRAM' in unit_label and any((v in unit_label for v in ['M**3', 'M^3', 'M3'])), 'Unverified concentration units: ' + unit_label)
            outputs = [v for v in data['output_files'] if Path(v['path']).suffix.lower() == '.csv']
            require(len(outputs) == 1, 'Expected one daily CSV')
            record = outputs[0]
            path = Path(record['path'])
            require(path.parent == folder and (not path.is_symlink()), 'Unexpected CSV path')
            stat = path.stat()
            require(stat.st_size == record['size'] and stat.st_mtime_ns == record['mtime_ns'], 'CSV changed: ' + str(path))
            require(list(pd.read_csv(path, nrows=0).columns) == ['utm_easting', 'utm_northing'] + columns, 'CSV calendar/schema mismatch')
            seen = set()
            rows = 0
            for chunk in pd.read_csv(path, chunksize=1000):
                xy = chunk[['utm_easting', 'utm_northing']].round(2)
                tuples = list(map(tuple, xy.to_numpy()))
                require(len(set(tuples)) == len(tuples) and (not seen.intersection(tuples)), 'Repeated receptors in ' + str(path))
                seen.update(tuples)
                indices = ENV['index'].get_indexer(pd.MultiIndex.from_frame(xy))
                require((indices >= 0).all(), 'Unmapped receptor coordinates in ' + str(path))
                values = chunk[columns].to_numpy(dtype=float)
                require(np.isfinite(values).all() and (values >= 0).all(), 'Invalid/missing concentrations in ' + str(path))
                accumulated[indices] += values
                rows += len(chunk)
            require(rows == data['unique_receptors'], 'Receptor count mismatch')
            included += 1
            points_processed += rows
            met_review += bool(data['model'].get('missing_met_hours_reported', 0))
            if source in ENV.get('source_review',{}).get(workflow,set()):
                source_review.append(source)
        grid_sums = np.zeros((len(ENV['gids']), days), dtype=float)
        positions = ENV.get('positions', np.arange(len(ENV['codes'])))
        for start in range(0, len(ENV['codes']), 5000):
            codes = ENV['codes'][start:start + 5000]
            valid = codes >= 0
            receptor_positions = positions[start:start + 5000][valid]
            np.add.at(grid_sums, codes[valid], np.asarray(accumulated[receptor_positions]))
        counts = ENV['counts']
        grid_means = np.full_like(grid_sums, np.nan)
        np.divide(grid_sums, counts[:, None], out=grid_means, where=counts[:, None] > 0)
        destination.mkdir(parents=True, exist_ok=True)
        daily_tmp = temporary / 'daily.csv.gz'
        annual_tmp = temporary / 'annual.csv'
        pd.DataFrame(grid_means, index=pd.Index(ENV['gids'], name='grid_ID'), columns=columns).to_csv(daily_tmp, compression='gzip')
        annual_values = grid_means.mean(axis=1)
        pd.DataFrame(dict(grid_ID=ENV['gids'], annual_mean_conc=annual_values, reference_receptors=counts, valid_days=np.where(counts > 0, days, 0))).to_csv(annual_tmp, index=False)
        shutil.copyfile(daily_tmp, daily)
        shutil.copyfile(annual_tmp, annual)
        result = dict(**config, state='completed_with_review_flags' if source_review or not (counts > 0).all() else 'completed', unresolved_source_assignment_review=source_review, positive_inventory_sources=int(positive.sum()), included_sources=included, boundary_excluded_source_ids=excluded, zero_emission_sources=int((~positive).sum()), modeled_receptor_rows=points_processed, sources_with_documented_missing_met_hours=met_review, reference_receptors=len(ENV['index']), output_grids=len(counts), grids_without_reference_receptors=[g for g, n in zip(ENV['gids'], counts) if n == 0], units='micrograms/m3, inherited from AERMOD CONC output; verify model unit labels before publication', aggregation='sum source concentrations once per unique receptor, then equal-weight mean selected receptors for each grid; receptor reuse across grids is allowed; annual arithmetic mean of calendar-day concentrations', assumptions='Unmodeled source contributions beyond each run receptor footprint are zero under existing spatial truncation. Missing expected in-boundary positive-emission cases block the group. Missing reference grids stay NaN, not zero.', scope='category-specific only; no UOGD cross-category total or source fractions', inventory=str(source_inventory), reference_mapping=str(REFERENCE_MAPPING), mapping=str(MAPPING), daily=str(daily), annual=str(annual))
        audit.write_text(json.dumps(result, indent=2))
        return result
    finally:
        if db:
            db.close()
        if accumulated is not None:
            accumulated.flush()
            del accumulated
        shutil.rmtree(temporary)

def safe_group(task):
    try:
        return compile_group(task)
    except Exception as exc:
        return dict(workflow=task[0], year=task[1], pollutant=task[2], state='blocked', error=str(exc))

def main():
    require(MAX_WORKERS >= 1 and MAX_WORKERS <= int(os.environ.get('SLURM_CPUS_PER_TASK', str(MAX_WORKERS))), 'Check CPU allocation')
    for p in [MAPPING, REFERENCE_MAPPING, GRID, BOUNDARY]:
        require(p.is_file(), 'Required compilation input missing: ' + str(p))
    env = load_environment()
    print('Reference receptors:', len(env['index']), 'output grids:', len(env['gids']), flush=True)
    del env
    OUTPUT.mkdir(parents=True, exist_ok=True)
    None
    summary = Path(_configured_path('compilation_summary'))
    tasks = [(w, y, p) for w in WORKFLOWS for y in YEARS for p in POLLUTANTS]
    if RETRY_SUMMARY is not None:
        with Path(RETRY_SUMMARY).open() as previous:
            records = [json.loads(line) for line in previous]
        tasks = sorted({(r['workflow'], int(r['year']), r['pollutant']) for r in records if r['state'] == 'blocked'})
        require(tasks, 'No blocked groups in retry summary')
    counts = Counter()
    with summary.open('x') as log, ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=initialize) as pool:
        futures = [pool.submit(safe_group, task) for task in tasks]
        for future in as_completed(futures):
            result = future.result()
            counts[result['state']] += 1
            log.write(json.dumps(result) + '\n')
            log.flush()
            print(result['workflow'], result['year'], result['pollutant'], result['state'], result.get('error', ''), flush=True)
    print('RESULT:', dict(counts), 'summary:', summary, flush=True)
    require(not counts['blocked'], 'Some groups blocked; inspect summary. Completed compiled groups are retained.')
if __name__ == '__main__':
    main()
