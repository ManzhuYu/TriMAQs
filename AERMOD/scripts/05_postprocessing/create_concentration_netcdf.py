"""Combine UOGD source categories into daily and annual NetCDF concentrations."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting, relative_metadata
from pathlib import Path
from datetime import datetime
import calendar, json, hashlib, gc
import numpy as np
import pandas as pd
from netCDF4 import Dataset
BASE = Path(_configured_path('create_concentration_netcdf.base_1'))
COMPILED = BASE / 'outputs/compiled'
SUMMARY = Path(_configured_path('compilation_summary'))
GRID = BASE / 'inputs/supporting/grids/study_grid_1km.csv'
SPEC = BASE / 'inputs/supporting/ml_grid/grid_spec_and_mask.npz'
OUTPUT = BASE / 'outputs/concentrations'
YEARS = _setting('create_concentration_netcdf.years')
POLLUTANTS = _setting('create_concentration_netcdf.pollutants')
SOURCES = ['Tier_1', 'Tier_2', 'UOGD']
WORKFLOWS = ['drilling_completion', 'production', 'compressors', 'other_tier2']
FACTORS = {'NO2': 24.45 / 46.01, 'SO2': 24.45 / 64.07, 'PM2.5': 1.0, 'VOC': 1.0}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def approved_absence(workflow, year, pollutant):
    return workflow == 'drilling_completion' and year < 2006 or (workflow == 'production' and (year < 2007 or pollutant == 'SO2')) or (workflow == 'compressors' and year < 2007)

def grid_environment():
    with np.load(SPEC, allow_pickle=False) as spec:
        data = {k: spec[k].copy() for k in ['x', 'y', 'lat', 'lon', 'study_mask']}
        data['crs'] = str(spec['crs'])
    grid = pd.read_csv(GRID)
    ids = grid.grid_ID.to_numpy(dtype=np.int64)
    ny, nx = (len(data['y']), len(data['x']))
    mask = data['study_mask']
    require(mask.dtype == bool and mask.shape == (ny, nx), 'Invalid study mask')
    require(len(np.unique(ids)) == len(ids) and (ids >= 0).all() and (ids < ny * nx).all(), 'Invalid grid IDs')
    require(set(ids) == set(np.flatnonzero(~mask.ravel())), 'Grid IDs differ from unmasked raster cells')
    require(np.allclose(data['lat'].ravel()[ids], grid.lat, rtol=0, atol=5e-06) and np.allclose(data['lon'].ravel()[ids], grid.lon, rtol=0, atol=5e-06), 'Grid ID/raster coordinate alignment differs')
    data.update(ids=ids, index=pd.Index(ids), ny=ny, nx=nx)
    return data

def add_category(target, record, env, year, pollutant):
    """Read each CSV once, require exact grid/calendar coverage and annual consistency."""
    workflow = record['workflow']
    state = record['state']
    if state == 'inventory_absent':
        require(approved_absence(workflow, year, pollutant), 'Unapproved absent inventory: ' + str(record))
        return dict(workflow=workflow, state='zero_by_scope_assumption', inventory=record.get('path'))
    require(state in {'completed', 'completed_with_review_flags'}, 'Category not successfully compiled: ' + str(record))
    daily = Path(record['daily'])
    annual = Path(record['annual'])
    require(daily.parent == COMPILED / workflow / pollutant and annual.parent == daily.parent, 'Unexpected category path')
    cols = [f'day{i}_conc' for i in range(1, target.shape[1] + 1)]
    require(list(pd.read_csv(daily, nrows=0).columns) == ['grid_ID'] + cols, 'Daily schema/calendar mismatch: ' + str(daily))
    annual_frame = pd.read_csv(annual).set_index('grid_ID')
    require(not annual_frame.index.duplicated().any() and set(annual_frame.index) == set(env['ids']), 'Annual grid mismatch')
    seen = np.zeros(len(env['ids']), dtype=bool)
    largest_difference = 0.0
    for chunk in pd.read_csv(daily, chunksize=2000):
        positions = env['index'].get_indexer(chunk.grid_ID)
        require((positions >= 0).all() and len(np.unique(positions)) == len(positions) and (not seen[positions].any()), 'Unknown/duplicate daily grid')
        values = chunk[cols].to_numpy(dtype=np.float64)
        require(np.isfinite(values).all() and (values >= 0).all(), 'Nonfinite/negative category concentrations')
        expected = annual_frame.loc[chunk.grid_ID, 'annual_mean_conc'].to_numpy(dtype=float)
        calculated = values.mean(axis=1)
        require(np.isfinite(expected).all() and np.allclose(calculated, expected, rtol=1e-09, atol=1e-12), 'Daily/annual category mismatch')
        largest_difference = max(largest_difference, float(np.max(abs(calculated - expected))))
        target[positions] += values
        seen[positions] = True
    require(seen.all(), 'Missing grids in daily category CSV')
    return dict(workflow=workflow, state=state, daily=str(daily), annual=str(annual), maximum_daily_annual_difference_ug_m3=largest_difference, source_assignment_review=record.get('unresolved_source_assignment_review', []))

def initialize_file(path, env, year, provenance):
    ds = Dataset(path, 'w', format='NETCDF4')
    days = 366 if calendar.isleap(year) else 365
    for name, size in [('time', days), ('source', 3), ('y', env['ny']), ('x', env['nx'])]:
        ds.createDimension(name, size)
    ds.title = f'Daily and annual UOGD AERMOD concentrations ({year})'
    ds.crs = env['crs']
    ds.grid_coordinate_note = 'x and y are raster indices; lat/lon provide geographic coordinates. CRS label inherited; full projection parameters are not supplied.'
    ds.source_definition = 'Tier_1=drilling/completion+production; Tier_2=compressors+other UOGD; UOGD=Tier_1+Tier_2. Do not sum all three source entries.'
    ds.aggregation = 'Equal-weight mean of four current AERMAP receptors per grid after summing source contributions; annual arithmetic mean of daily concentrations.'
    ds.gas_conversion = 'NO2 and SO2: ppb=(ug/m3)*24.45/molecular_weight; reference molar volume at 25 degC and 1 atm, matching prior export.'
    ds.processing_summary = relative_metadata(str(SUMMARY))
    ds.summary_sha256 = provenance
    ds.created = datetime.now().isoformat()
    ds.assumptions = 'Known pre-activity inventory gaps and production SO2 excluded from inventory contribute zero. Existing model footprint truncation retained. Some production source assignments remain flagged in input audits.'
    time = ds.createVariable('time', 'i4', ('time',))
    time[:] = np.arange(days)
    time.units = f'days since {year}-01-01 00:00:00'
    time.calendar = 'proleptic_gregorian'
    source = ds.createVariable('source', str, ('source',))
    source[:] = np.array(SOURCES, dtype=object)
    for name in ['x', 'y']:
        v = ds.createVariable(name, 'i4', (name,))
        v[:] = env[name]
        v.long_name = 'raster ' + ('column' if name == 'x' else 'row') + ' index'
    for name, unit in [('lat', 'degrees_north'), ('lon', 'degrees_east')]:
        v = ds.createVariable(name, 'f4', ('y', 'x'), zlib=True, complevel=4)
        v[:] = env[name]
        v.units = unit
        v.standard_name = 'latitude' if name == 'lat' else 'longitude'
    mask = ds.createVariable('study_mask', 'i1', ('y', 'x'), zlib=True, complevel=4)
    mask[:] = env['study_mask']
    mask.long_name = '1 outside study area; 0 inside'
    gid = ds.createVariable('grid_ID', 'i4', ('y', 'x'), fill_value=-1, zlib=True, complevel=4)
    raster = np.full((env['ny'], env['nx']), -1, dtype=np.int32)
    raster.ravel()[env['ids']] = env['ids']
    gid[:] = raster
    return ds

def write_pollutant(ds, env, pol, tier1, tier2):
    name = pol.replace('.', '')
    days = tier1.shape[1]
    shape = (env['ny'], env['nx'])
    factor = FACTORS[pol]
    unit = 'ppb' if pol in {'NO2', 'SO2'} else 'ug/m3'
    daily = ds.createVariable(name, 'f4', ('time', 'source', 'y', 'x'), fill_value=np.nan, zlib=True, complevel=4, chunksizes=(1, 1, *shape))
    annual = ds.createVariable(name + '_annual', 'f4', ('source', 'y', 'x'), fill_value=np.nan, zlib=True, complevel=4, chunksizes=(1, *shape))
    for v in [daily, annual]:
        v.units = unit
        v.coordinates = 'lat lon'
    daily.long_name = 'Daily ' + name + ' concentration'
    annual.long_name = 'Annual mean ' + name + ' concentration'
    daily.cell_methods = 'time: mean'
    annual.annual_mean_definition = 'Arithmetic mean of all daily means in this calendar year'
    for source in range(3):
        mean = np.zeros(len(env['ids']), dtype=np.float64)
        for day in range(days):
            values = (tier1[:, day] if source == 0 else tier2[:, day] if source == 1 else tier1[:, day] + tier2[:, day]) * factor
            values = values.astype(np.float32)
            require(np.isfinite(values).all(), 'float32 overflow')
            raster = np.full(shape, np.nan, dtype=np.float32)
            raster.ravel()[env['ids']] = values
            daily[day, source] = raster
            mean += values
        raster = np.full(shape, np.nan, dtype=np.float32)
        raster.ravel()[env['ids']] = mean / days
        annual[source] = raster

def verify_file(path, env, year):
    days = 366 if calendar.isleap(year) else 365
    maximum_error = 0.0
    with Dataset(path) as ds:
        require(len(ds.dimensions['time']) == days, 'Wrong saved calendar')
        require(list(ds.variables['source'][:]) == SOURCES, 'Wrong saved source labels')
        for pol in POLLUTANTS:
            name = pol.replace('.', '')
            daily = ds.variables[name]
            mean = np.zeros((3, len(env['ids'])), dtype=float)
            for day in range(days):
                slices = np.ma.filled(daily[day], np.nan)
                require(np.isnan(slices[:, env['study_mask']]).all(), 'Outside mask is not NaN')
                values = slices.reshape(3, -1)[:, env['ids']]
                require(np.isfinite(values).all() and (values >= 0).all(), 'Missing/invalid study concentrations')
                require(np.allclose(values[2], values[0] + values[1], rtol=2e-06, atol=1e-10), 'Saved UOGD daily total differs from tier sum')
                mean += values
            annual = np.ma.filled(ds.variables[name + '_annual'][:], np.nan)
            require(np.isnan(annual[:, env['study_mask']]).all(), 'Annual outside mask is not NaN')
            annual = annual.reshape(3, -1)[:, env['ids']]
            require(np.allclose(annual, mean / days, rtol=2e-06, atol=1e-10), 'Saved annual mean differs from saved daily means')
            require(np.allclose(annual[2], annual[0] + annual[1], rtol=2e-06, atol=1e-10), 'Saved annual tier sum mismatch')
            maximum_error = max(maximum_error, float(np.max(abs(annual - mean / days))))
    return maximum_error

def main():
    env = grid_environment()
    raw = SUMMARY.read_bytes()
    rows = [json.loads(line) for line in raw.decode().splitlines() if line.strip()]
    records = {(r['workflow'], int(r['year']), r['pollutant']): r for r in rows}
    require(len(records) == len(rows), 'Duplicate compilation summary identities')
    for year in YEARS:
        for pol in POLLUTANTS:
            for workflow in WORKFLOWS:
                key = (workflow, year, pol)
                require(key in records, 'Missing summary group: ' + str(key))
                r = records[key]
                require(r['state'] in {'completed', 'completed_with_review_flags'} or (r['state'] == 'inventory_absent' and approved_absence(*key)), 'Unapproved summary state: ' + str(r))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    log = OUTPUT / ('EXPORT_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '.jsonl')
    failed = 0
    sha = hashlib.sha256(raw).hexdigest()
    with log.open('x') as handle:
        for year in YEARS:
            final = OUTPUT / f'aermod_daily_conc_{year}.nc'
            temporary = OUTPUT / f'aermod_daily_conc_{year}.partial.nc'
            try:
                require(not final.exists() and (not temporary.exists()), 'Output/partial exists; inspect before rerunning: ' + str(final))
                audit = []
                days = 366 if calendar.isleap(year) else 365
                with initialize_file(temporary, env, year, sha) as ds:
                    for pol in POLLUTANTS:
                        print(year, pol, 'combining', flush=True)
                        first = np.zeros((len(env['ids']), days), dtype=np.float64)
                        second = np.zeros_like(first)
                        for workflow, target in zip(WORKFLOWS, [first, first, second, second]):
                            audit.append(dict(pollutant=pol, **add_category(target, records[workflow, year, pol], env, year, pol)))
                        write_pollutant(ds, env, pol, first, second)
                        del first, second
                        gc.collect()
                    ds.category_audit = json.dumps(relative_metadata(audit))
                print(year, 'checking saved NetCDF', flush=True)
                difference = verify_file(temporary, env, year)
                temporary.rename(final)
                result = dict(year=year, state='completed', output=str(final), size=final.stat().st_size, maximum_annual_readback_difference=difference, categories=audit)
            except Exception as exc:
                failed += 1
                result = dict(year=year, state='blocked', error=str(exc), partial=str(temporary))
            handle.write(json.dumps(relative_metadata(result)) + '\n')
            handle.flush()
            print(year, result['state'], result.get('error', ''), flush=True)
    print('Share:', log, flush=True)
    require(failed == 0, 'Some years blocked; see export log. Existing/partial files retained.')
if __name__ == '__main__':
    main()
