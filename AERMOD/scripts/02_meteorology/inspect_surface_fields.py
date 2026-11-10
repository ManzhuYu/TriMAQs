"""Inspect surface meteorological fields for missing or invalid values."""
from __future__ import annotations
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
from collections import OrderedDict
from contextlib import contextmanager
import csv
from datetime import datetime, timedelta
import hashlib
import json
import logging
from pathlib import Path
import platform
import shutil
import tempfile
import zipfile
import numpy as np
years = _setting('inspect_surface_fields.years')
process_kinds = ['surface', 'upper']
resume_completed = True
grid_csv_file = _configured_path('inspect_surface_fields.grid_csv_file_1')
surface_base_input_dir = _configured_path('inspect_surface_fields.surface_base_input_dir_2')
land_input_dir = _configured_path('inspect_surface_fields.land_input_dir_3')
cloud_input_dir = _configured_path('inspect_surface_fields.cloud_input_dir_4')
boundary_input_dir = _configured_path('inspect_surface_fields.boundary_input_dir_5')
surface_output_dir = _configured_path('inspect_surface_fields.surface_output_dir_6')
upper_output_dir = _configured_path('inspect_surface_fields.upper_output_dir_7')
upper_input_dir = _configured_path('inspect_surface_fields.upper_input_dir_8')
extra_upper_files = [boundary_input_dir + '/pressure_level_data_4hrs_01Jan24.nc']
summary_file = surface_output_dir + '/REBUILD_SUMMARY.json'
extra_land_files = [boundary_input_dir + '/Land_hourly_4hrs_01Jan24/data_0.nc', boundary_input_dir + '/Land_hourly_4hrs_01Jan24/data_1.nc']
extra_cloud_files = [boundary_input_dir + '/total_cloud_cover_4hrs_01Jan24.nc']
negative_roundoff_tolerance = {'tp': 2e-06, 'ssrd': 32.0, 'strd': 32.0}
grid_ids = None
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
ACCUMULATED = {'tp', 'ssrd', 'strd'}
SURFACE_UNITS = {'t2m': 'K', 'd2m': 'K', 'u10': 'm s-1', 'v10': 'm s-1', 'sp': 'Pa', 'tp': 'm', 'tcc': '1', 'ssrd': 'J m-2', 'strd': 'J m-2'}
UPPER_UNITS = {'t': 'K', 'u': 'm s-1', 'v': 'm s-1', 'z': 'm2 s-2', 'r': '%'}
ALIASES = {'K': {'k', 'kelvin'}, 'm s-1': {'ms-1', 'ms**-1', 'm/s'}, 'Pa': {'pa', 'pascal', 'pascals'}, 'm': {'m', 'metres', 'meters'}, '1': {'1', '(0-1)', 'fraction'}, 'J m-2': {'jm-2', 'jm**-2', 'j/m2'}, 'm2 s-2': {'m2s-2', 'm**2s**-2', 'm2/s2'}, '%': {'%', 'percent'}}

def require(condition, message):
    if not condition:
        raise ValueError(message)

def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=str) + '\n', encoding='utf-8')

def script_provenance():
    """Notebook cells have no __file__; hashing must never stop processing."""
    name = globals().get('__file__')
    if not name:
        return {'reader_script_sha256': None, 'execution_mode': 'notebook_or_interactive', 'provenance_note': 'No __file__ in this session; save the notebook/cell with the run records.'}
    try:
        return {'reader_script_sha256': hashlib.sha256(Path(name).read_bytes()).hexdigest(), 'execution_mode': 'python_file', 'reader_script_path': str(Path(name).resolve())}
    except OSError as exc:
        return {'reader_script_sha256': None, 'execution_mode': 'file_path_unavailable', 'provenance_note': str(exc)}

def attributes(obj):
    return {n: str(obj.getncattr(n)) for n in obj.ncattrs()}

def load_grids(path, selected=None):
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        row['grid_ID'] = str(int(float(row['grid_ID'])))
    require(len({r['grid_ID'] for r in rows}) == len(rows), 'Duplicate grid IDs')
    if selected is not None:
        ids = {str(int(g)) for g in selected}
        rows = [r for r in rows if r['grid_ID'] in ids]
        require({r['grid_ID'] for r in rows} == ids, 'Some selected grid IDs were not found')
    require(bool(rows), 'No grids selected')
    return rows

def source_members(path):
    path = Path(path)
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            members = sorted((i.filename for i in z.infolist() if Path(i.filename).suffix.lower() in ('.nc', '.nc4', '.cdf')))
        require(bool(members), f'No NetCDF members in {path}')
        return [(str(path), member) for member in members]
    return [(str(path), None)]

@contextmanager
def open_source(spec):
    import netCDF4 as nc
    path, member = spec
    if member is None:
        with nc.Dataset(path) as ds:
            yield ds
    else:
        with zipfile.ZipFile(path) as z:
            if z.getinfo(member).file_size <= 128 * 1024 * 1024:
                payload = z.read(member)
                with nc.Dataset('raw_memory.nc', memory=payload) as ds:
                    yield ds
            else:
                with tempfile.TemporaryDirectory(prefix='era5_raw_') as temporary:
                    temp_path = Path(temporary) / Path(member).name
                    with z.open(member) as src, temp_path.open('wb') as dst:
                        shutil.copyfileobj(src, dst, length=1024 * 1024)
                    with nc.Dataset(temp_path) as ds:
                        yield ds

class DatasetPool:
    """Small serial LRU; no shared netCDF handles across threads/processes."""

    def __init__(self, maximum=4):
        self.maximum, self.items = (maximum, OrderedDict())

    def get(self, spec):
        if spec in self.items:
            self.items.move_to_end(spec)
            return self.items[spec][1]
        while len(self.items) >= self.maximum:
            _, (context, _) = self.items.popitem(last=False)
            context.__exit__(None, None, None)
        context = open_source(spec)
        ds = context.__enter__()
        self.items[spec] = (context, ds)
        return ds

    def close(self):
        while self.items:
            _, (context, _) = self.items.popitem()
            context.__exit__(None, None, None)

def regular_weights(lat, lon, target_lat, target_lon):
    """Positive bilinear weights on regular 1D coordinates; no extrapolation."""
    lat, lon = (np.asarray(lat, dtype=float), (np.asarray(lon, dtype=float) + 180) % 360 - 180)
    require(lat.ndim == lon.ndim == 1 and len(lat) >= 2 and (len(lon) >= 2), 'Expected 1D regular lat/lon coordinates')
    require(np.isfinite(lat).all() and np.isfinite(lon).all(), 'Missing spatial coordinates')
    ai, bi = (np.argsort(lat), np.argsort(lon))
    a, b = (lat[ai], lon[bi])
    require((np.diff(a) > 0).all() and (np.diff(b) > 0).all(), 'Duplicate spatial coordinates or unsupported longitude seam')
    x, y = (np.asarray(target_lat), (np.asarray(target_lon) + 180) % 360 - 180)
    require(((x >= a[0]) & (x <= a[-1]) & (y >= b[0]) & (y <= b[-1])).all(), 'Grid centre outside raw domain; no nearest-neighbour extrapolation is applied')
    ia = np.clip(np.searchsorted(a, x, side='right') - 1, 0, len(a) - 2)
    ib = np.clip(np.searchsorted(b, y, side='right') - 1, 0, len(b) - 2)
    wa, wb = ((x - a[ia]) / (a[ia + 1] - a[ia]), (y - b[ib]) / (b[ib + 1] - b[ib]))
    corners = [(ai[ia], bi[ib]), (ai[ia + 1], bi[ib]), (ai[ia], bi[ib + 1]), (ai[ia + 1], bi[ib + 1])]
    weights = [(1 - wa) * (1 - wb), wa * (1 - wb), (1 - wa) * wb, wa * wb]
    return (corners, weights)

def apply_weights(field, weights):
    corners, fractions = weights
    values = np.stack([field[..., a, b] for a, b in corners])
    w = np.asarray(fractions)
    while w.ndim < values.ndim:
        w = np.expand_dims(w, axis=1)
    valid = (np.isfinite(values) | (w <= 1e-12)).all(axis=0)
    answer = (np.where(np.isfinite(values), values, 0) * w).sum(axis=0)
    return np.where(valid, answer, np.nan)

def hourly_increment(now, previous, when, tolerance):
    """Traditional ERA5-Land: 01 UTC is step 1; 00 UTC is preceding step 24."""
    increment = np.asarray(now, dtype=float).copy() if when.hour == 1 else np.asarray(now, dtype=float) - np.asarray(previous, dtype=float) if previous is not None else np.full_like(now, np.nan, dtype=float)
    tiny = (increment < 0) & (increment >= -tolerance)
    large = increment < -tolerance
    increment[tiny] = 0
    increment[large] = np.nan
    return (increment, int(tiny.sum()), int(large.sum()))

class Catalog:
    """Attach every field to its original decoded timestamp before any sorting."""

    def __init__(self, units, grids, upper=False):
        self.units, self.upper = (units, upper)
        self.series = {n: {} for n in units}
        self.providers, self.sources, self.weights = ([], [], {})
        self.levels = None
        self.lat = np.array([float(g['center_lat']) for g in grids])
        self.lon = np.array([float(g['center_lon']) for g in grids])
        self.pool = DatasetPool()
        self.cache = OrderedDict()

    def add(self, path, allowed, product):
        import netCDF4 as nc
        path = Path(path)
        require(path.is_file(), f'Missing configured raw input: {path}')
        for spec in source_members(path):
            with open_source(spec) as ds:
                names = [n for n in allowed if n in ds.variables]
                if not names:
                    continue
                tn = next((n for n in ('valid_time', 'time') if n in ds.variables), None)
                require(tn is not None, f'{spec}: no time coordinate')
                tv = ds.variables[tn]
                require(len(tv.dimensions) == 1, f'{spec}: multidimensional forecast time is unsupported')
                require(getattr(tv, 'calendar', 'standard') in ('standard', 'gregorian', 'proleptic_gregorian'), 'Unsupported calendar')
                decoded = nc.num2date(tv[:], tv.units, calendar=getattr(tv, 'calendar', 'standard'))
                times = [datetime(t.year, t.month, t.day, t.hour, t.minute, t.second, t.microsecond) for t in decoded]
                require(times and all((t.minute == t.second == t.microsecond == 0 for t in times)), 'Empty or off-hour raw time coordinates')
                lan = next((n for n in ('latitude', 'lat') if n in ds.variables), None)
                lon = next((n for n in ('longitude', 'lon') if n in ds.variables), None)
                require(lan is not None and lon is not None, f'{spec}: missing coordinates')
                require(len(ds.variables[lan].dimensions) == len(ds.variables[lon].dimensions) == 1, f'{spec}: requires one-dimensional latitude/longitude coordinates')
                ld, od = (ds.variables[lan].dimensions[0], ds.variables[lon].dimensions[0])
                require(ld != od, 'Raw input must be a regular lat/lon grid, not paired locations')
                levn, level_order = (None, None)
                if self.upper:
                    levn = next((n for n in ('pressure_level', 'level') if n in ds.variables), None)
                    require(levn is not None, f'{spec}: no pressure levels')
                    lv = ds.variables[levn]
                    lu = str(getattr(lv, 'units', '')).lower().replace(' ', '')
                    require(lu in ('hpa', 'pa', 'millibars', 'millibar', 'mb'), f'{spec}: ambiguous pressure-level units {lu!r}')
                    levels = np.asarray(lv[:], dtype=float) / (100 if lu == 'pa' else 1)
                    require(np.isfinite(levels).all() and len(set(levels)) == len(levels) and (levels > 0).all(), 'Invalid pressure levels')
                    level_order = np.argsort(-levels)
                    levels = levels[level_order]
                    require(self.levels is None or np.array_equal(levels, self.levels), 'Pressure levels differ across raw inputs')
                    self.levels = levels
                    levn = lv.dimensions[0]
                provider = {'spec': spec, 'time_dim': tv.dimensions[0], 'lat_dim': ld, 'lon_dim': od, 'level_dim': levn, 'level_order': level_order}
                provider_id = len(self.providers)
                self.providers.append(provider)
                self.weights[provider_id] = regular_weights(ds.variables[lan][:], ds.variables[lon][:], self.lat, self.lon)
                checked = {}
                for name in names:
                    var = ds.variables[name]
                    unit = str(getattr(var, 'units', '')).lower().replace(' ', '')
                    require(unit in ALIASES[self.units[name]], f'{spec}: {name} units {unit!r}; expected {self.units[name]}. No raw-unit guessing.')
                    expected_dims = {provider['time_dim'], ld, od}
                    if self.upper:
                        expected_dims.add(levn)
                    require(expected_dims.issubset(set(var.dimensions)), f'{spec}: missing expected dimensions for {name}: {var.dimensions}')
                    extra = set(var.dimensions) - expected_dims
                    require(all((len(ds.dimensions[d]) == 1 for d in extra)), f'{spec}: non-singleton ensemble/forecast/expver axis in {name}; select/merge explicitly upstream')
                    if name in ACCUMULATED:
                        require(product in ('traditional_era5_land', 'declared_cumulative_01UTC'), 'Accumulated fields require an explicitly configured cumulative source')
                        require(str(getattr(var, 'GRIB_stepType', 'accum')) == 'accum', f'{spec}: {name} is not declared accumulated')
                    checked[name] = {'dimensions': list(var.dimensions), 'attributes': attributes(var)}
                    for i, dt in enumerate(times):
                        self.series[name].setdefault(dt, []).append((provider_id, name, i))
                stat = path.stat()
                self.sources.append({'path': str(path), 'archive_member': spec[1], 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns, 'product': product, 'first_stored': str(times[0]), 'last_stored': str(times[-1]), 'minimum': str(min(times)), 'maximum': str(max(times)), 'records': len(times), 'backward_steps': sum((b < a for a, b in zip(times, times[1:]))), 'global_attributes': attributes(ds), 'variables': checked})
                logging.info('Indexed %s%s: %s times; variables=%s', path.name, f'::{spec[1]}' if spec[1] else '', len(times), names)

    def raw_points(self, reference):
        pid, name, i = reference
        provider = self.providers[pid]
        ds = self.pool.get(provider['spec'])
        var = ds.variables[name]
        retained, selections = ([], [])
        for dim in var.dimensions:
            if dim == provider['time_dim']:
                selections.append(i)
            elif dim in (provider['lat_dim'], provider['lon_dim'], provider['level_dim']):
                selections.append(slice(None))
                retained.append(dim)
            else:
                selections.append(0)
        field = np.ma.asarray(var[tuple(selections)]).astype(float).filled(np.nan)
        desired = ([provider['level_dim']] if self.upper else []) + [provider['lat_dim'], provider['lon_dim']]
        field = np.transpose(field, [retained.index(d) for d in desired])
        if self.upper:
            field = field[provider['level_order']]
        return apply_weights(field, self.weights[pid])

    def points(self, name, when):
        key = (name, when)
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        references = self.series[name].get(when, [])
        shape = (len(self.levels), len(self.lat)) if self.upper else (len(self.lat),)
        result = np.full(shape, np.nan)
        for j, reference in enumerate(references):
            candidate = self.raw_points(reference)
            if j:
                require(np.allclose(result, candidate, rtol=1e-06, atol=1e-06, equal_nan=True), f"Conflicting duplicate {name} at {when}: sources={[self.providers[r[0]]['spec'] for r in references]}; no silent overwrite")
            else:
                result = candidate
        self.cache[key] = result
        while len(self.cache) > 32:
            self.cache.popitem(last=False)
        return result

    def close(self):
        self.pool.close()

def required_times(year, upper=False):
    start, end = (datetime(year, 1, 1), datetime(year + 1, 1, 1))
    n = int((end - start).total_seconds() / 3600)
    if upper:
        return [start + timedelta(hours=h) for h in range(5, n + 5) if (start + timedelta(hours=h)).hour in (0, 6, 12, 18)]
    return [start + timedelta(hours=h) for h in range(6, n + 6)]

def write_processed(catalog, grids, year, destination, tolerances=None):
    import netCDF4 as nc
    destination = Path(destination)
    partial = destination.with_name(destination.name + '.partial')
    manifest_file = destination.with_suffix('.manifest.json')
    require(not destination.exists() and (not partial.exists()) and (not manifest_file.exists()), f'Output/partial/manifest already exists for {destination}; choose a new output_dir')
    destination.parent.mkdir(parents=True, exist_ok=True)
    start, end = (datetime(year, 1, 1), datetime(year + 1, 1, 1) + timedelta(hours=5))
    times = sorted({t for records in catalog.series.values() for t in records if start <= t <= end})
    require(bool(times), 'No source timestamps in requested year/boundary window')
    required = set(required_times(year, catalog.upper))
    missing_times = sorted(required - set(times))
    units = catalog.units
    required_variables = [n for n in units if catalog.upper or n != 'strd']
    missing_by_variable = {n: [str(t) for t in sorted(required - set(catalog.series[n]))] for n in units}
    for name, absent in missing_by_variable.items():
        if absent:
            logging.warning('%s: %s required raw timestamps missing; first 12: %s', name, len(absent), absent[:12])
    problems = {n: {'required_times_with_nonfinite_values': [], 'nonfinite_grid_values_required': 0, 'roundoff_negative_increments_clipped': 0, 'large_negative_increments': 0, 'missing_previous_accumulation_times': []} for n in units}
    manifest = {'state': 'started', 'year': year, 'output': str(destination), 'time_basis': 'UTC', 'local_calendar': 'EST UTC-5, no DST', 'spatial_interpolation': 'positive bilinear weights; no extrapolation or masked-corner filling', 'units': units, 'source_files': catalog.sources, 'grid_csv_order': [g['grid_ID'] for g in grids], 'source_time_policy': 'Each value attached to decoded raw UTC time before sorting; duplicates checked', 'accumulation_policy': 'tp/ssrd/strd declared cumulative with 01 UTC forecast-day reset; folder name alone is not product identification', 'negative_roundoff_tolerance': tolerances or {}, 'python_version': platform.python_version(), 'netCDF4_version': nc.__version__, 'numpy_version': np.__version__, 'records': len(times), 'first_UTC': str(times[0]), 'last_UTC': str(times[-1]), 'required_timestamp_count': len(required), 'missing_required_UTC_times': [str(t) for t in missing_times]}
    manifest['missing_required_raw_times_by_variable'] = missing_by_variable
    manifest.update(script_provenance())
    manifest['aermet_required_variables'] = required_variables
    manifest['supporting_variables'] = [n for n in units if n not in required_variables]
    dump(manifest_file, manifest)
    try:
        with nc.Dataset(partial, 'w', format='NETCDF4') as out:
            out.createDimension('time', len(times))
            out.createDimension('location', len(grids))
            if catalog.upper:
                require(catalog.levels is not None, 'No upper-air levels found')
                out.createDimension('level', len(catalog.levels))
                lv = out.createVariable('level', 'f4', ('level',))
                lv[:] = catalog.levels
                lv.units, lv.positive, lv.standard_name = ('hPa', 'down', 'air_pressure')
            else:
                out.createDimension('bnds', 2)
            tv = out.createVariable('time', 'f8', ('time',))
            tv.units, tv.calendar = (f'hours since {year}-01-01 00:00:00', 'proleptic_gregorian')
            tv.standard_name, tv.timezone = ('time', 'UTC')
            encoded = nc.date2num(times, tv.units, calendar=tv.calendar)
            tv[:] = encoded
            if not catalog.upper:
                bounds = out.createVariable('time_bnds', 'f8', ('time', 'bnds'))
                bounds[:] = np.column_stack([encoded - 1, encoded])
                tv.bounds = 'time_bnds'
            location = out.createVariable('location', 'i4', ('location',))
            location[:] = np.arange(len(grids))
            ids = out.createVariable('grid_ID', 'i8', ('location',))
            ids[:] = [int(g['grid_ID']) for g in grids]
            for name, values, unit in (('latitude', catalog.lat, 'degrees_north'), ('longitude', catalog.lon, 'degrees_east')):
                var = out.createVariable(name, 'f8', ('location',))
                var[:], var.units, var.standard_name = (values, unit, name)
            variables = {}
            for name, unit in units.items():
                dims = ('level', 'time', 'location') if catalog.upper else ('time', 'location')
                chunks = (1, min(24, len(times)), min(256, len(grids))) if catalog.upper else (min(24, len(times)), min(256, len(grids)))
                var = out.createVariable(name, 'f4', dims, zlib=True, complevel=2, chunksizes=chunks, fill_value=np.float32(np.nan))
                var.set_var_chunk_cache(size=16 * 1024 * 1024, nelems=1009, preemption=0.75)
                var.units, var.coordinates = (unit, 'latitude longitude')
                var.cell_methods = 'time: sum' if name in ACCUMULATED else 'time: point'
                if name in ACCUMULATED:
                    var.processing = 'Hourly increment ending at time, deaccumulated in UTC with ERA5-Land step-1 reset at 01 UTC'
                variables[name] = var
            out.Conventions, out.time_basis = ('CF-1.8', 'UTC')
            out.history = 'Rebuilt from raw UTC records; timestamps and values sorted together; never trimmed by record count.'
            out.spatial_interpolation = manifest['spatial_interpolation']
            out.product_status = 'in_progress'
            for i, dt in enumerate(times):
                for name in units:
                    values = catalog.points(name, dt)
                    if name in ACCUMULATED:
                        previous_time = dt - timedelta(hours=1)
                        previous = None if dt.hour == 1 or previous_time not in catalog.series[name] else catalog.points(name, previous_time)
                        if dt.hour != 1 and previous is None:
                            problems[name]['missing_previous_accumulation_times'].append(str(dt))
                        values, tiny, large = hourly_increment(values, previous, dt, tolerances[name])
                        problems[name]['roundoff_negative_increments_clipped'] += tiny
                        problems[name]['large_negative_increments'] += large
                    if dt in required:
                        missing = int((~np.isfinite(values)).sum())
                        if missing:
                            problems[name]['required_times_with_nonfinite_values'].append(str(dt))
                            problems[name]['nonfinite_grid_values_required'] += missing
                    if catalog.upper:
                        variables[name][:, i, :] = values
                    else:
                        variables[name][i, :] = values
                if i % 168 == 0 or i == len(times) - 1:
                    logging.info('%s: wrote %s/%s UTC hours (%s)', destination.name, i + 1, len(times), dt)
                    out.sync()
            ready = not missing_times and all((not problems[n]['required_times_with_nonfinite_values'] for n in required_variables))
            out.product_status = 'ready_for_AERMET_validation_checks' if ready else 'incomplete_required_meteorology'
            out.aermet_required_input_complete = str(ready).lower()
            out.aermet_required_calendar = 'four UTC synoptic soundings per EST day' if catalog.upper else 'EST hour-ending annual calendar'
            out.setncattr('missing_required_UTC_timestamps', json.dumps([str(t) for t in missing_times]))
        partial.rename(destination)
        manifest.update(state='completed', required_input_complete=ready, variable_diagnostics=problems, warning='Completion means preprocessing finished, not scientific/model validation.')
        dump(manifest_file, manifest)
        logging.info('SAVED %s | required AERMET input complete=%s', destination, ready)
        if not ready:
            logging.warning('NOT READY FOR ANNUAL AERMET: inspect %s; missing timestamps=%s; first 12=%s', manifest_file, len(missing_times), [str(t) for t in missing_times[:12]])
        return manifest
    except Exception as exc:
        manifest.update(state='failed', error=str(exc), variable_diagnostics=problems)
        dump(manifest_file, manifest)
        raise
    finally:
        catalog.close()

def surface_sources(year):
    """Monthly files mandatory; previous December/next January are optional auxiliaries."""
    sources = []
    base, land, cloud = (Path(surface_base_input_dir), Path(land_input_dir), Path(cloud_input_dir))
    for month in range(1, 13):
        sources.extend([(base / f'{year}{month:02d}.nc', 'base'), (land / f'{year}_{month:02d}_era5_land.nc', 'land'), (cloud / f'tcc_{year}_{month:02d}.nc', 'cloud')])
    for path, kind in [(base / f'{year - 1}12.nc', 'base'), (base / f'{year + 1}01.nc', 'base'), (land / f'{year - 1}_12_era5_land.nc', 'land'), (land / f'{year + 1}_01_era5_land.nc', 'land'), (cloud / f'tcc_{year + 1}_01.nc', 'cloud')]:
        if path.is_file():
            sources.append((path, kind))
    if year == 2023:
        for name, kind in [(p, 'boundary') for p in extra_land_files] + [(p, 'cloud') for p in extra_cloud_files]:
            if Path(name).is_file():
                sources.append((Path(name), kind))
            else:
                logging.warning('Optional 2024 boundary file unavailable: %s', name)
    return list(dict.fromkeys(sources))

def index_surface_sources(catalog, sources):
    policies = {'base': (['t2m', 'd2m', 'u10', 'v10', 'sp', 'tp'], 'declared_cumulative_01UTC'), 'land': (['ssrd', 'strd'], 'traditional_era5_land'), 'cloud': (['tcc'], 'ERA5_cloud'), 'boundary': ([n for n in SURFACE_UNITS if n != 'tcc'], 'traditional_era5_land')}
    for path, kind in sources:
        names, convention = policies[kind]
        catalog.add(path, names, convention)

def check_main_year_variables(catalog, year):
    """Stop before writing if a field exists only in a short boundary file."""
    start, end = (datetime(year, 1, 1), datetime(year + 1, 1, 1))
    absent = [n for n in SURFACE_UNITS if n != 'strd' and (not any((start <= t < end for t in catalog.series[n])))]
    require(not absent, f'No main-year raw records for {absent}. Check the monthly source mapping; boundary-only values cannot replace the year.')

def finite_json(values):
    return [float(v) if np.isfinite(v) else None for v in np.asarray(values).reshape(-1)]

def describe_corners(catalog, name, when, indices):
    reports = []
    for pid, variable, i in catalog.series[name].get(when, []):
        provider = catalog.providers[pid]
        ds = catalog.pool.get(provider['spec'])
        var = ds.variables[variable]
        selections = []
        retained = []
        for dim in var.dimensions:
            if dim == provider['time_dim']:
                selections.append(i)
            elif dim in (provider['lat_dim'], provider['lon_dim']):
                selections.append(slice(None))
                retained.append(dim)
            else:
                selections.append(0)
        field = np.ma.asarray(var[tuple(selections)]).astype(float).filled(np.nan)
        field = np.transpose(field, [retained.index(provider['lat_dim']), retained.index(provider['lon_dim'])])
        corners, weights = catalog.weights[pid]
        points = []
        for j in indices:
            points.append({'location_index': int(j), 'corners': [{'row': int(a[j]), 'column': int(b[j]), 'weight': float(w[j]), 'raw_value': float(field[a[j], b[j]]) if np.isfinite(field[a[j], b[j]]) else None} for (a, b), w in zip(corners, weights)]})
        reports.append({'source': provider['spec'], 'variable_attributes': attributes(var), 'points': points})
    return reports

def main():
    grids = load_grids(grid_csv_file, grid_ids)
    reports = []
    manifests = sorted(Path(surface_output_dir).glob('*_surface_variables_yearly.manifest.json'))
    require(manifests, 'No surface manifests found; check surface_output_dir')
    for path in manifests:
        manifest = json.loads(path.read_text())
        if manifest.get('required_input_complete'):
            continue
        result = {'manifest': str(path), 'year': manifest['year'], 'missing_raw_times_by_variable': manifest.get('missing_required_raw_times_by_variable'), 'variable_diagnostics': manifest.get('variable_diagnostics'), 'events': []}
        for name, diagnostic in manifest['variable_diagnostics'].items():
            times = [datetime.fromisoformat(t) for t in diagnostic['required_times_with_nonfinite_values']]
            if not times:
                continue
            needed = set(times) | {t - timedelta(hours=1) for t in times}
            units = {name: SURFACE_UNITS[name]}
            catalog = Catalog(units, grids)
            try:
                added = set()
                for info in manifest['source_files']:
                    if name not in info['variables']:
                        continue
                    low = datetime.fromisoformat(info['minimum'])
                    high = datetime.fromisoformat(info['maximum'])
                    if not any((low <= t <= high for t in needed)):
                        continue
                    key = (info['path'], info['product'])
                    if key not in added:
                        catalog.add(Path(info['path']), [name], info['product'])
                        added.add(key)
                for when in times:
                    now = catalog.points(name, when)
                    previous_time = when - timedelta(hours=1)
                    previous = catalog.points(name, previous_time) if previous_time in catalog.series[name] else None
                    if name in ACCUMULATED:
                        increment, tiny, large = hourly_increment(now, previous, when, negative_roundoff_tolerance[name])
                    else:
                        increment, tiny, large = (now, 0, 0)
                    bad = np.flatnonzero(~np.isfinite(increment))
                    selected = bad[:20]
                    event = {'variable': name, 'time_UTC': str(when), 'affected_grids': len(bad), 'large_negative_increments': large, 'small_negative_increments_clipped': tiny, 'sample_grid_IDs': [grids[int(i)]['grid_ID'] for i in selected], 'current_interpolated_values': finite_json(now[selected]), 'previous_interpolated_values': finite_json(previous[selected]) if previous is not None else None, 'current_raw_corners': describe_corners(catalog, name, when, selected), 'previous_raw_corners': describe_corners(catalog, name, previous_time, selected)}
                    result['events'].append(event)
                    logging.info('Checked %s %s: affected grids=%s', name, when, len(bad))
            except Exception as exc:
                result['events'].append({'variable': name, 'diagnostic_error': str(exc)})
                logging.exception('Diagnostic failed for %s %s', manifest['year'], name)
            finally:
                catalog.close()
        reports.append(result)
    destination = Path(_configured_path('surface_inspection_report'))
    destination.parent.mkdir(parents=True,exist_ok=True)
    require(not destination.exists(), 'Report already exists')
    dump(destination, reports)
    print('Surface inspection report: ' + str(destination))
    print('Existing NetCDF, raw files and manifests were not modified.')
if __name__ == '__main__':
    main()
