"""Run both AERMET stages for selected years and meteorological grids."""
from __future__ import annotations
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
import logging
import calendar
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback
import numpy as np
grid_csv_file = _configured_path('run_aermet.grid_csv_file_1')
surface_nc_file_folder = _configured_path('run_aermet.surface_nc_file_folder_2')
upperair_nc_file_folder = _configured_path('run_aermet.upperair_nc_file_folder_3')
aersurface_output_dir = _configured_path('run_aermet.aersurface_output_dir_4')
aermap_output_dir = _configured_path('run_aermet.aermap_output_dir_5')
aermet_executable = _configured_path('run_aermet.aermet_executable_6')
output_dir = _configured_path('run_aermet.output_dir_7')
years = _setting('run_aermet.years')
period_start = None
period_end = None
grid_ids = [33980]
max_workers = _setting('run_aermet.max_workers')
resume_completed = True
recheck_previous_script_outputs = True
archive_pressure_rounding_failures = True
archive_unfinished_cases_from_stopped_job = True
archive_confirmed_upperair_elevation_failures = True
archive_orphan_folders_from_stopped_job = True
inspect_only = False
prepare_only = False
surface_time_basis = 'UTC'
upper_time_basis = 'UTC'
local_utc_offset_hours = -5
radiation_mode = 'hourly_energy'
precipitation_mode = 'hourly_depth'
surface_additional_files = []
upper_additional_files = []
elevation_column = None
expected_aermet_version = '24142'
max_location_distance_km = 2.0
minimum_profile_depth_m = 3000
maximum_noncalm_invalid_hours = 0
units_if_missing = {}
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def build_configuration():
    """Collect the settings above for worker processes; no external config file."""
    return {'grid_csv': grid_csv_file, 'period_start': period_start, 'period_end': period_end, 'surface_pattern': str(Path(surface_nc_file_folder) / '{year}_surface_variables_yearly.nc'), 'upper_pattern': str(Path(upperair_nc_file_folder) / '{year}_upper_air_final.nc'), 'aersurface_pattern': str(Path(aersurface_output_dir) / '{year}' / 'Grid_{grid_id}' / 'aersurface_out_{grid_id}_{year}.sfc'), 'terrain_pattern': str(Path(aermap_output_dir) / 'Grid_{grid_id}' / 'source_output_{grid_id}_grid.DAT'), 'output_root': output_dir, 'aermet_executable': aermet_executable, 'surface_time_basis': surface_time_basis, 'upper_time_basis': upper_time_basis, 'local_utc_offset_hours': local_utc_offset_hours, 'radiation_mode': radiation_mode, 'precipitation_mode': precipitation_mode, 'surface_additional_files': surface_additional_files, 'upper_additional_files': upper_additional_files, 'elevation_column': elevation_column, 'expected_aermet_version': expected_aermet_version, 'max_location_distance_km': max_location_distance_km, 'minimum_profile_depth_m': minimum_profile_depth_m, 'maximum_noncalm_invalid_hours': maximum_noncalm_invalid_hours, 'units_if_missing': units_if_missing}
SURFACE = ('t2m', 'd2m', 'u10', 'v10', 'sp', 'tp', 'tcc', 'ssrd')
UPPER = ('u', 'v', 't', 'z', 'r')
MONTHS = ('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC')
UNIT_OPTIONS = {'K': {'k', 'kelvin'}, 'm/s': {'ms-1', 'ms**-1', 'm/s', 'm.s-1'}, 'Pa': {'pa', 'pascal', 'pascals'}, 'm': {'m', 'metre', 'meter', 'metres', 'meters'}, '1': {'1', '(0-1)', 'fraction'}, 'J/m2': {'jm-2', 'jm**-2', 'j/m2'}, 'W/m2': {'wm-2', 'wm**-2', 'w/m2'}, '%': {'%', 'percent'}, 'm2/s2': {'m2s-2', 'm**2s**-2', 'm^2s^-2', 'm2/s2'}, 'hPa': {'hpa', 'millibars', 'millibar', 'mb'}}

def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, default=str) + '\n', encoding='utf-8')

def require(condition, message):
    if not condition:
        raise ValueError(message)

def calculate_wind_direction(u, v):
    return math.degrees(math.atan2(-u, -v)) % 360 if u or v else 0.0

def rh_from_dewpoint(t, d):
    t, d = (t - 273.15, d - 273.15)
    return min(100.0, max(0.0, 100 * math.exp(17.625 * d / (243.04 + d) - 17.625 * t / (243.04 + t))))

def dewpoint(t, rh):
    require(0 < rh <= 100, f'Invalid upper-air RH: {rh}; zero RH cannot define a dew point')
    a = math.log(rh / 100) + 17.625 * (t - 273.15) / (243.04 + t - 273.15)
    return 243.04 * a / (17.625 - a) + 273.15

def humidity_to_fsl_dewpoint(temp_k, rh, audit=None, context=None):
    """Unsupported humidity produces a missing dew point, not fabricated RH.

    FSL version 2 accepts 99999 as missing. Values >100 are outside this validation's
    conversion range; that is not a claim that all supersaturation is unphysical.
    The existing approximate Magnus conversion for supported RH is retained.
    ERA5 mixed-phase humidity needs separate quantitative thermodynamic review.
    """
    reason = None
    if not math.isfinite(rh):
        reason = 'nonfinite_RH'
    elif rh <= 0:
        reason = 'nonpositive_RH_no_finite_dewpoint'
    elif rh > 100:
        reason = 'RH_above_supported_conversion_range'
    if reason:
        event = dict(context or {})
        event.update(raw_RH_percent=rh, temperature_K=temp_k, action='FSL_dewpoint_missing_99999', reason=reason)
        if audit is not None:
            audit.append(event)
        return 99999
    return round((dewpoint(temp_k, rh) - 273.15) * 10)

def period_bounds(year, cfg=None):
    cfg = cfg or {}
    a, b = (cfg.get('period_start'), cfg.get('period_end'))
    require(bool(a) == bool(b), 'Set both period_start and period_end, or both None')
    start = datetime.fromisoformat(a) if a else datetime(year, 1, 1)
    last = datetime.fromisoformat(b) if b else datetime(year, 12, 31)
    require(start.year == last.year == year and start <= last, 'Validation dates must be ordered within the selected year')
    return (start, last + timedelta(days=1))

def hourly_times(year, cfg=None):
    """Local hour-ending 1..24 for the selected complete calendar days."""
    start, end = period_bounds(year, cfg)
    count = int((end - start).total_seconds() / 3600)
    return [start + timedelta(hours=h) for h in range(1, count + 1)]

def script_hash():
    name = globals().get('__file__')
    return hashlib.sha256(Path(name).read_bytes()).hexdigest() if name and Path(name).is_file() else None

def to_local(dt, basis, offset):
    require(basis in ('UTC', 'local_standard'), 'Set each time_basis to UTC or local_standard after inspecting the processed files')
    return dt + timedelta(hours=offset) if basis == 'UTC' else dt

def named_slice(var, location, upper=False):
    """Do not infer axis order from shape. Index each variable by named axes."""
    dims = list(var.dimensions)
    require(dims.count('location') == 1 and dims.count('time') == 1, f"{getattr(var, 'name', 'variable')}: expected location/time axes; got {dims}")
    retained = [d for d in dims if d != 'location']
    wanted = ['level', 'time'] if upper else ['time']
    require(sorted(retained) == sorted(wanted), f'Unsupported dimensions {dims}; expected {wanted} plus location')
    index = [location if d == 'location' else slice(None) for d in dims]
    data = np.ma.asarray(var[tuple(index)]).astype(float).filled(np.nan)
    return np.transpose(data, [retained.index(d) for d in wanted])

def check_unit(var, expected, fallback, provenance):
    actual = str(getattr(var, 'units', '')).strip()
    unit = actual or fallback
    require(bool(unit), f'{var.name}: no units attribute; supply an explicitly verified units_if_missing entry')
    normal = unit.lower().replace(' ', '')
    require(normal in UNIT_OPTIONS[expected], f'{var.name}: units {unit!r}, expected {expected}; convert upstream or extend the reader explicitly')
    provenance[var.name] = {'units': unit, 'source': 'NetCDF attribute' if actual else 'configured fallback'}

def file_info(path):
    p = Path(path)
    stat = p.stat()
    return {'path': str(p.resolve()), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}

def input_paths(cfg, kind, year):
    paths = [Path(cfg[kind + '_pattern'].format(year=year))]
    paths += [Path(p.format(year=year, previous_year=year - 1, next_year=year + 1)) for p in cfg.get(kind + '_additional_files', [])]
    for p in paths:
        require(p.is_file(), f'Missing input: {p}')
    return list(dict.fromkeys(paths))

def nearest_location(ds, lat, lon, tolerance):
    latvar = next((ds.variables[n] for n in ('latitude', 'lat') if n in ds.variables), None)
    lonvar = next((ds.variables[n] for n in ('longitude', 'lon') if n in ds.variables), None)
    require(latvar is not None and lonvar is not None, 'No latitude/longitude coordinate variables')
    require(latvar.dimensions == ('location',) and lonvar.dimensions == ('location',), 'This reader requires paired one-dimensional location coordinates; regular lat/lon grids need an explicit adapter')
    la = np.ma.asarray(latvar[:]).astype(float).filled(np.nan)
    lo = (np.ma.asarray(lonvar[:]).astype(float).filled(np.nan) + 180) % 360 - 180
    require(np.isfinite(la).all() and np.isfinite(lo).all(), 'Masked/nonfinite location coordinates')
    dy = np.radians(la - lat)
    dx = np.radians((lo - lon + 180) % 360 - 180)
    a = np.sin(dy / 2) ** 2 + math.cos(math.radians(lat)) * np.cos(np.radians(la)) * np.sin(dx / 2) ** 2
    km = 6371 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
    idx = int(np.argmin(km))
    require(km[idx] <= tolerance, f'Nearest NetCDF location is {km[idx]:.2f} km away; limit={tolerance}')
    return (idx, {'index': idx, 'latitude': float(la[idx]), 'longitude': float(lo[idx]), 'distance_km': float(km[idx])})

def read_series(cfg, kind, year, lat, lon):
    import netCDF4 as nc
    upper = kind == 'upper'
    fields = UPPER if upper else SURFACE
    expected = dict(zip(SURFACE, ('K', 'K', 'm/s', 'm/s', 'Pa', 'm', '1', 'J/m2')))
    if cfg.get('radiation_mode') == 'hourly_mean_flux':
        expected['ssrd'] = 'W/m2'
    expected.update({'u': 'm/s', 'v': 'm/s', 't': 'K', 'z': 'm2/s2', 'r': '%'})
    records, provenance, reference_levels = ({}, [], None)
    window_start = datetime(year, 1, 1) - timedelta(days=2)
    window_end = datetime(year + 1, 1, 1) + timedelta(days=2)
    for path in input_paths(cfg, kind, year):
        with nc.Dataset(str(path)) as ds:
            idx, location = nearest_location(ds, lat, lon, cfg['max_location_distance_km'])
            require(all((n in ds.variables for n in fields)), f'{path}: required variables: {fields}')
            tv = ds.variables['time']
            require(tv.dimensions == ('time',), 'Time variable must use the time dimension')
            require(getattr(tv, 'calendar', 'standard') in ('standard', 'gregorian', 'proleptic_gregorian'), 'Unsupported calendar')
            dates = nc.num2date(tv[:], tv.units, calendar=getattr(tv, 'calendar', 'standard'))
            require(len(dates) > 0, f'{path}: empty time coordinate')
            times = [to_local(datetime(d.year, d.month, d.day, d.hour, d.minute, d.second), cfg[kind + '_time_basis'], cfg['local_utc_offset_hours']) for d in dates]
            require(all((t.minute == t.second == 0 for t in times)), 'Non-hourly timestamp alignment')
            require(all((b > a for a, b in zip(times, times[1:]))), f'{path}: duplicate or unsorted time coordinates')
            units = {}
            arrays = {}
            for n in fields:
                check_unit(ds.variables[n], expected[n], cfg['units_if_missing'].get(n), units)
                arrays[n] = named_slice(ds.variables[n], idx, upper)
                require(arrays[n].shape[-1] == len(times), f'{path}: {n} time axis mismatch')
            if upper:
                lv = ds.variables['level']
                check_unit(lv, 'hPa', cfg['units_if_missing'].get('level'), units)
                levels = np.asarray(lv[:], dtype=float)
                require(np.isfinite(levels).all() and len(set(levels)) == len(levels) and (levels > 0).all(), 'Invalid pressure levels')
                order = np.argsort(-levels)
                levels = levels[order]
                arrays = {n: a[order, :] for n, a in arrays.items()}
                require(reference_levels is None or np.array_equal(reference_levels, levels), 'Pressure levels change across input files')
                reference_levels = levels
            for j, dt in enumerate(times):
                if not window_start <= dt <= window_end:
                    continue
                row = {n: a[:, j].copy() if upper else float(a[j]) for n, a in arrays.items()}
                if dt in records:
                    require(all((np.allclose(row[n], records[dt][n], rtol=0, atol=1e-06, equal_nan=True) for n in fields)), f'Conflicting overlapping {kind} data at {dt}; correct the input file list')
                else:
                    records[dt] = row
            provenance.append({**file_info(path), 'location': location, 'units': units, 'first_local_time': str(times[0]), 'last_local_time': str(times[-1]), 'variables': {n: {'dimensions': list(ds.variables[n].dimensions), 'shape': list(ds.variables[n].shape)} for n in fields}})
    return (records, reference_levels, provenance)

def require_coverage(records, expected, label):
    missing = [d for d in expected if d not in records]
    require(not missing, f'{label}: {len(missing)} missing timestamps; first 12: {missing[:12]}. Supply genuine boundary/source data; no zero filling is performed.')

def surface_values(row, cfg):
    require(all((math.isfinite(row[n]) for n in SURFACE)), 'Missing/nonfinite surface data')
    require(180 < row['t2m'] < 340 and 180 < row['d2m'] < 340, 'Invalid Kelvin surface temperature')
    require(row['d2m'] <= row['t2m'] + 0.5, 'Dew point exceeds temperature by >0.5 K')
    require(50000 < row['sp'] < 110000 and 0 <= row['tcc'] <= 1, 'Invalid pressure/cloud fraction')
    require(-1e-10 <= row['tp'] <= 0.254, 'Invalid hourly precipitation depth (metres)')
    require(cfg['precipitation_mode'] == 'hourly_depth', 'Precipitation must be confirmed as hourly depth; deaccumulate upstream first')
    mode = cfg['radiation_mode']
    require(mode in ('hourly_energy', 'hourly_mean_flux'), 'Set radiation_mode to hourly_energy or hourly_mean_flux; cumulative forecasts require upstream deaccumulation')
    radiation = row['ssrd'] / 3600 if mode == 'hourly_energy' else row['ssrd']
    require(-0.1 <= radiation <= 1400, 'Invalid hourly solar radiation; check accumulation convention')
    speed = math.hypot(row['u10'], row['v10'])
    require(speed < 100, 'Invalid 10 m wind speed')
    return (calculate_wind_direction(row['u10'], row['v10']), speed, row['t2m'] - 273.15, row['sp'] / 10, rh_from_dewpoint(row['t2m'], row['d2m']), max(0, radiation), row['tcc'] * 10, max(0, row['tp']) * 100000)

def create_onsite_met(path, records, year, cfg):
    times = hourly_times(year, cfg)
    require_coverage(records, times, 'Onsite local hour-ending coverage')
    with Path(path).open('w', encoding='ascii', newline='\n') as f:
        for dt in times:
            day = dt - timedelta(days=1) if dt.hour == 0 else dt
            hour = dt.hour or 24
            vals = surface_values(records[dt], cfg)
            f.write(f'{day.year} {day.month} {day.day} {hour}\n')
            f.write(' '.join((f'{v:.5f}' for v in vals[:5])) + '\n')
            f.write(' '.join((f'{v:.5f}' for v in vals[5:])) + '\n')

def profile_levels(row, levels, surface, elevation, min_depth, humidity_audit=None, context=None, geometry_audit=None):
    p0 = surface['sp'] / 100
    output = [[9, round(p0 * 10), round(elevation), round((surface['t2m'] - 273.15) * 10), round((surface['d2m'] - 273.15) * 10), round(calculate_wind_direction(surface['u10'], surface['v10'])), round(math.hypot(surface['u10'], surface['v10']) * 10)]]
    previous_h, previous_p = (round(elevation), round(p0 * 10))
    for i, p in enumerate(levels):
        raw_z = float(row['z'][i])
        height_m = raw_z / 9.80665 if math.isfinite(raw_z) else float('nan')
        height = round(height_m) if math.isfinite(height_m) else None
        reason = None
        if p >= p0:
            reason = 'pressure_at_or_below_model_surface'
        elif round(float(p) * 10) == round(p0 * 10):
            reason = 'pressure_indistinguishable_from_surface_at_FSL_precision'
        elif height is not None and (height_m <= elevation or height <= round(elevation)):
            reason = 'at_or_below_selected_terrain_reference'
        if reason:
            event = dict(context or {})
            event.update(pressure_hPa=float(p), surface_pressure_hPa=p0, selected_surface_elevation_m=elevation, raw_geopotential_m2_s2=raw_z, geopotential_height_m=height_m, rounded_height_m=height, action='omit_level_from_FSL', reason=reason)
            if geometry_audit is not None:
                geometry_audit.append(event)
            continue
        require(all((math.isfinite(float(row[n][i])) for n in UPPER if n != 'r')), f'Missing upper-air value at {p} hPa; time={context}')
        require(height > previous_h, f'Nonincreasing ABOVE-REFERENCE heights: level={p} hPa, height={height_m:.3f} m, previous={previous_h} m, terrain={elevation:.3f} m, surface_pressure={p0:.3f} hPa; time={context}')
        pressure = round(float(p) * 10)
        require(pressure < previous_p, f'Nondecreasing rounded profile pressure: level={p} hPa, encoded={pressure}, previous_encoded={previous_p}, surface={p0} hPa; time={context}')
        t, rh, u, v = (float(row[n][i]) for n in ('t', 'r', 'u', 'v'))
        require(150 < t < 340 and math.hypot(u, v) < 200, 'Invalid upper-air temperature/wind')
        humidity_context = dict(context or {}, pressure_hPa=float(p))
        dewpoint_code = humidity_to_fsl_dewpoint(t, rh, humidity_audit, humidity_context)
        output.append([4, pressure, height, round((t - 273.15) * 10), dewpoint_code, round(calculate_wind_direction(u, v)), round(math.hypot(u, v) * 10)])
        previous_h, previous_p = (height, pressure)
    require(len(output) >= 4 and output[-1][2] - elevation >= min_depth, f'Insufficient above-ground profile: {len(output)} levels, top {output[-1][2] - elevation:.0f} m AGL')
    return output

def convert_to_fsl(path, upper, levels, surface, year, lat, lon, elevation, cfg):
    offset = cfg['local_utc_offset_hours']
    local_start, local_end = period_bounds(year, cfg)
    utc_start = local_start - timedelta(hours=offset)
    utc_end = local_end - timedelta(hours=offset)
    expected = []
    dt = utc_start.replace(hour=0)
    while dt < utc_end:
        if dt >= utc_start and dt.hour in (0, 6, 12, 18):
            expected.append(dt + timedelta(hours=offset))
        dt += timedelta(hours=1)
    require_coverage(upper, expected, 'Upper-air synoptic coverage')
    require_coverage(surface, expected, 'Surface values for FSL soundings')
    tops, humidity_audit = ([], [])
    pressure_level_count = 0
    retained_counts = []
    geometry_audit = []
    geometry_path = Path(path).with_name('upperair_geometry_audit.csv')
    audit_path = Path(path).with_name('upperair_humidity_audit.csv')
    columns = ['time_UTC', 'time_EST', 'pressure_hPa', 'temperature_K', 'raw_RH_percent', 'action', 'reason']
    geometry_columns = ['time_UTC', 'time_EST', 'pressure_hPa', 'surface_pressure_hPa', 'selected_surface_elevation_m', 'raw_geopotential_m2_s2', 'geopotential_height_m', 'rounded_height_m', 'action', 'reason']
    with Path(path).open('w', encoding='ascii', newline='\n') as f, audit_path.open('w', encoding='utf-8', newline='') as audit_file, geometry_path.open('w', encoding='utf-8', newline='') as geometry_file:
        audit_writer = csv.DictWriter(audit_file, fieldnames=columns)
        audit_writer.writeheader()
        geometry_writer = csv.DictWriter(geometry_file, fieldnames=geometry_columns)
        geometry_writer.writeheader()
        for local in expected:
            utc = local - timedelta(hours=offset)
            surface_values(surface[local], cfg)
            audit_start = len(humidity_audit)
            geometry_start = len(geometry_audit)
            try:
                data = profile_levels(upper[local], levels, surface[local], elevation, cfg['minimum_profile_depth_m'], humidity_audit, {'time_UTC': str(utc), 'time_EST': str(local)}, geometry_audit)
            finally:
                audit_writer.writerows(humidity_audit[audit_start:])
                audit_file.flush()
                geometry_writer.writerows(geometry_audit[geometry_start:])
                geometry_file.flush()
            pressure_level_count += len(data) - 1
            retained_counts.append(len(data) - 1)
            f.write(f'254 {utc.hour} {utc.day} {MONTHS[utc.month - 1]} {utc.year}\n')
            f.write(f"1 99999 99999 {abs(lat):.3f}{('N' if lat >= 0 else 'S')} {abs(lon):.3f}{('E' if lon >= 0 else 'W')}\n")
            f.write(f'2 99999 99999 99999 {len(data) + 4}\n')
            f.write(f'{3:7d}{99999:7d}{99999:7d}{99999:7d}'.ljust(47) + 'ms\n')
            for row in data:
                line = ''.join((f'{v:7d}' for v in row))
                require(len(line) <= 60, 'FSL record exceeds AERMET input width')
                f.write(line + '\n')
            tops.append(data[-1][2] - elevation)
    pressure_excluded = sum((r['reason'] == 'pressure_at_or_below_model_surface' for r in geometry_audit))
    terrain_excluded = len(geometry_audit) - pressure_excluded
    if terrain_excluded:
        logging.warning('Excluded %s pressure levels at/below selected terrain; review %s before accepting the validation.', terrain_excluded, geometry_path)
    affected = len({r['time_UTC'] for r in humidity_audit})
    if humidity_audit:
        logging.warning('Upper-air dew point marked missing at %s/%s pressure levels in %s/%s soundings. Review %s', len(humidity_audit), pressure_level_count, affected, len(expected), audit_path)
    return {'soundings': len(expected), 'minimum_top_agl_m': min(tops), 'maximum_top_agl_m': max(tops), 'pressure_level_count': pressure_level_count, 'missing_dewpoint_levels': len(humidity_audit), 'soundings_with_missing_dewpoint': affected, 'humidity_audit_csv': str(audit_path), 'humidity_policy': 'Unsupported RH produces missing FSL dew point; raw data unchanged; no clipping/floor', 'geometry_audit_csv': str(geometry_path), 'pressure_excluded_levels': pressure_excluded, 'terrain_reference_excluded_levels': terrain_excluded, 'minimum_retained_pressure_levels_per_sounding': min(retained_counts), 'geometry_policy': 'Require pressure below surface pressure AND height above selected site terrain; no height shifting or invented levels'}

def generate_aermet_stage1_inp(folder, grid_id, year, lat, lon, elevation, cfg):
    """Write the original-style stage-1 filename with corrected units/heights/clocks."""
    xy = f"{abs(lon):.6f}{('E' if lon >= 0 else 'W')} {abs(lat):.6f}{('N' if lat >= 0 else 'S')}"
    first, end = period_bounds(year, cfg)
    date = f'{first:%Y/%m/%d} TO {end - timedelta(days=1):%Y/%m/%d}'
    adjust = -cfg['local_utc_offset_hours']
    stage1 = f'JOB\n  MESSAGES AERMET_MSG.MES\n  REPORT AERMET_REPORT.RPT\nONSITE\n  DATA onsite_{grid_id}_{year}.met\n  QAOUT ONSITE_QA.OQA\n  XDATES {date}\n  LOCATION 99999 {xy} 0 {elevation:.2f}\n  OSHEIGHTS 2 10\n  READ 1 OSYR OSMO OSDY OSHR\n  FORMAT 1 FREE\n  READ 2 WD02 WS02 TT01 PRES RH01\n  FORMAT 2 FREE\n  READ 3 INSO TSKC PRCP\n  FORMAT 3 FREE\n  RANGE TT -90 <= 65 -9999\n  RANGE WS 0 <= 100 -9999\n  RANGE WD 0 <= 360 -9999\n  RANGE PRCP 0 <= 25400 -9999\n  RANGE RH 0 <= 100 -9999\n  RANGE INSO 0 <= 1400 -9999\n  RANGE PRES 5000 <= 11000 -9999\n  RANGE TSKC 0 <= 10 -9999\n  THRESHOLD 0.1\nUPPERAIR\n  DATA upperair_{grid_id}_{year}.fsl FSL\n  XDATES {date}\n  QAOUT UPPER_QA.OQA\n  EXTRACT UPPER_QA.IQA\n  LOCATION 99999 {xy} {adjust} {elevation:.2f}\n'
    path = folder / f'aermet_stage1_{grid_id}_{year}.inp'
    path.write_text(stage1, encoding='ascii')
    return path

def generate_aermet_stage2_inp(folder, grid_id, year, cfg=None):
    """Write the original-style stage-2 filename using local QA/AERSURFACE files."""
    first, end = period_bounds(year, cfg)
    date = f'{first:%Y/%m/%d} TO {end - timedelta(days=1):%Y/%m/%d}'
    stage2 = f'JOB\n  REPORT stage2.rpt\n  MESSAGES stage2.mes\nUPPERAIR\n  QAOUT UPPER_QA.OQA\nONSITE\n  QAOUT ONSITE_QA.OQA\nMETPREP\n  XDATES {date}\n  OUTPUT output_stage2.SFC\n  PROFILE output_stage2.PFL\n  AERSURF aersurface_out_{grid_id}_{year}.sfc\n'
    path = folder / f'aermet_stage2_{grid_id}_{year}.inp'
    path.write_text(stage2, encoding='ascii')
    return path

def completion_upperair_location_headers(folder, elevation):
    """Completion only the Fortran F7.3 overflow token; retain the true elevation.

    AERMET reads the LOCATION line as whitespace-delimited fields on Stage 2.
    This works around the QA writer's limited internal character field without
    modifying the executable or its upper-air profile rows.
    """
    changes = []
    require((folder / 'UPPER_QA.OQA').is_file(), 'Missing Stage-1 upper QA output')
    for name in ('UPPER_QA.OQA', 'UPPER_QA.IQA'):
        path = folder / name
        if not path.is_file():
            continue
        lines = path.read_text(encoding='ascii').splitlines(keepends=True)
        require(len(lines) > 2 and lines[1].startswith('LOCATION '), f'Unexpected upper-air QA header: {path}')
        tokens = lines[1].split()
        require(len(tokens) == 6, f'Expected LOCATION station lat lon offset elevation in {path}: {tokens}')
        before = lines[1].rstrip('\r\n')
        if set(tokens[-1]) == {'*'}:
            tokens[-1] = f'{elevation:.3f}'
            lines[1] = ' '.join(tokens) + '\n'
            path.write_text(''.join(lines), encoding='ascii')
            changes.append({'file': name, 'before': before, 'after': lines[1].rstrip(), 'true_elevation_m': elevation, 'reason': 'Fortran F7.3 station-elevation header overflow'})
        else:
            require(abs(float(tokens[-1]) - elevation) <= 0.0011, f'Upper QA elevation differs from site reference: {path}')
    dump(folder / 'upperair_header_format.json', {'changes': changes, 'true_site_elevation_m': elevation, 'profile_rows_changed': False})

def run_aermet_stage(exe, folder, stage, expected_version, inp_file, report_file):
    result = subprocess.run([str(Path(exe).resolve()), str(inp_file.name)], cwd=folder, capture_output=True, text=True, errors='replace', check=False)
    (folder / f'{stage}.stdout.log').write_text(result.stdout, encoding='utf-8')
    (folder / f'{stage}.stderr.log').write_text(result.stderr, encoding='utf-8')
    report = folder / report_file
    text = report.read_text(errors='replace') if report.exists() else ''
    combined = result.stdout + '\n' + result.stderr + '\n' + text
    require(result.returncode == 0, f'{stage} exited {result.returncode}; see logs')
    require('AERMET FINISHED SUCCESSFULLY' in combined.upper() and 'UN-SUCCESSFULLY' not in combined.upper(), f'{stage}: no unequivocal successful completion; see report/logs')
    require(re.search('(?<!\\d)' + re.escape(str(expected_version)) + '(?!\\d)', combined), f'{stage}: expected AERMET {expected_version} was not identified in report/logs')
    for name in ('ONSITE_QA.OQA', 'UPPER_QA.OQA') if stage == 'stage1' else ('output_stage2.SFC', 'output_stage2.PFL'):
        require((folder / name).is_file() and (folder / name).stat().st_size > 0, f'{stage}: missing {name}')

def output_key(parts, sfc):
    yr, mo, day = map(int, parts[:3])
    yr = 2000 + yr if yr < 70 else 1900 + yr if yr < 100 else yr
    hour = int(parts[4 if sfc else 3])
    require(1 <= hour <= 24, f'Invalid output hour: {hour}')
    return datetime(yr, mo, day) + timedelta(hours=hour)

def audit_outputs(folder, year, max_missing, cfg=None):
    expected = set(hourly_times(year, cfg))
    seen, invalid, calm = (set(), [], 0)
    with (folder / 'output_stage2.SFC').open() as f:
        header = next(f)
        require('VERSION:' in header, 'Missing SFC header')
        for line in f:
            p = line.split()
            if not p:
                continue
            require(len(p) >= 22, 'Truncated SFC row')
            key = output_key(p, True)
            require(key not in seen, f'Duplicate SFC hour {key}')
            seen.add(key)
            us, wsstar, zic, zim, length, speed, direction, temp = (float(p[i]) for i in (6, 7, 9, 10, 11, 15, 16, 18))
            bad = not all((math.isfinite(x) for x in (us, zim, length, speed, direction, temp)))
            bad |= speed < 0 or speed >= 90 or direction < 0 or (direction > 360) or (temp <= 0) or (temp > 900)
            bad |= us < 0 or us >= 9 or zim < 0 or (zim > 90000) or (length < -99990)
            bad |= length < 0 and (not math.isfinite(zic) or not math.isfinite(wsstar) or zic < 0 or (zic > 90000) or (wsstar < 0))
            if speed == 0:
                calm += 1
            elif bad:
                invalid.append(str(key))
    require(seen == expected, f'SFC calendar mismatch: {len(expected - seen)} missing, {len(seen - expected)} unexpected hours')
    pfl = set()
    with (folder / 'output_stage2.PFL').open() as f:
        for line in f:
            if line.strip():
                pfl.add(output_key(line.split(), False))
    require(pfl == expected, f'PFL calendar mismatch: {len(expected - pfl)} missing, {len(pfl - expected)} unexpected hours')
    result = {'hours': len(seen), 'calm_hours': calm, 'noncalm_invalid_hours': len(invalid), 'invalid_times': invalid, 'note': 'Automated screening only; inspect_netcdf_inputs AERMET warning/message files and validation results before production.'}
    dump(folder / 'output_audit.json', result)
    require(len(invalid) <= max_missing, f'{len(invalid)} noncalm hours fail meteorological screening (limit {max_missing}); inspect_netcdf_inputs output_audit.json')
    return result

def elevation_for(cfg, grid):
    column = cfg.get('elevation_column')
    if column:
        return (float(grid[column]), f'grid CSV column {column}')
    path = Path(cfg['terrain_pattern'].format(grid_id=grid['grid_ID']))
    values = []
    for line in path.read_text().splitlines():
        p = line.split()
        if len(p) >= 7 and p[0:2] == ['SO', 'LOCATION'] and (p[2] == str(grid['grid_ID'])):
            values.append(float(p[6]))
    require(len(values) == 1, f"{path}: expected one SO LOCATION for grid {grid['grid_ID']}")
    require(math.isfinite(values[0]), 'Invalid site elevation')
    return (values[0], str(path))

def process_grid(cfg, grid, year, prepare_only):
    gid = grid['grid_ID']
    folder = Path(cfg['output_root']) / str(year) / f'Grid_{gid}'
    require(not folder.exists(), f'{folder} already exists. Use a new output_dir; existing results are never overwritten.')
    folder.mkdir(parents=True)
    onsite_met_file = folder / f'onsite_{gid}_{year}.met'
    upper_air_fsl_file = folder / f'upperair_{gid}_{year}.fsl'
    aersurface_file = folder / f'aersurface_out_{gid}_{year}.sfc'
    logging.info('Processing Grid %s, Year %s -> %s', gid, year, folder)
    status = {'grid_id': gid, 'year': year, 'state': 'started', 'configuration': cfg, 'script_sha256': script_hash(), 'period_start': cfg.get('period_start'), 'period_end': cfg.get('period_end'), 'validation_only': bool(cfg.get('period_start'))}
    try:
        lat, lon = (float(grid['center_lat']), float(grid['center_lon']))
        elev, elev_source = elevation_for(cfg, grid)
        status.update(latitude=lat, longitude=lon, elevation_m=elev, elevation_source=elev_source)
        surface, _, sprov = read_series(cfg, 'surface', year, lat, lon)
        upper, levels, uprov = read_series(cfg, 'upper', year, lat, lon)
        status['source_files'] = {'surface': sprov, 'upper': uprov}
        dump(folder / 'manifest.json', status)
        create_onsite_met(onsite_met_file, surface, year, cfg)
        status['upperair_summary'] = convert_to_fsl(upper_air_fsl_file, upper, levels, surface, year, lat, lon, elev, cfg)
        source = Path(cfg['aersurface_pattern'].format(year=year, grid_id=gid))
        shutil.copyfile(source, aersurface_file)
        status['aersurface'] = file_info(source)
        stage1_inp_file = generate_aermet_stage1_inp(folder, gid, year, lat, lon, elev, cfg)
        stage2_inp_file = generate_aermet_stage2_inp(folder, gid, year, cfg)
        status['state'] = 'prepared_only'
        dump(folder / 'manifest.json', status)
        if not prepare_only:
            status['executable'] = file_info(cfg['aermet_executable'])
            run_aermet_stage(cfg['aermet_executable'], folder, 'stage1', cfg['expected_aermet_version'], stage1_inp_file, 'AERMET_REPORT.RPT')
            completion_upperair_location_headers(folder, elev)
            run_aermet_stage(cfg['aermet_executable'], folder, 'stage2', cfg['expected_aermet_version'], stage2_inp_file, 'stage2.rpt')
            status['output_audit'] = audit_outputs(folder, year, cfg['maximum_noncalm_invalid_hours'], cfg)
            status['output_files'] = {n: file_info(folder / n) for n in ('output_stage2.SFC', 'output_stage2.PFL')}
            status['state'] = 'completed_pending_scientific_review'
            dump(folder / 'DONE.json', status)
        dump(folder / 'manifest.json', status)
        return {'grid_id': gid, 'year': year, 'state': status['state'], 'folder': str(folder)}
    except Exception as exc:
        status.update(state='failed', error=str(exc))
        dump(folder / 'manifest.json', status)
        (folder / 'failure.txt').write_text(traceback.format_exc(), encoding='utf-8')
        raise

def inspect_netcdf_inputs(cfg, year):
    import netCDF4 as nc
    for kind in ('surface', 'upper'):
        for path in input_paths(cfg, kind, year):
            with nc.Dataset(path) as ds:
                print(f'\n{kind}: {path}')
                for name, var in ds.variables.items():
                    print(name, var.dimensions, var.shape, 'units=', getattr(var, 'units', '<missing>'))
                tv = ds.variables['time']
                dates = nc.num2date(tv[[0, len(tv) - 1]], tv.units, calendar=getattr(tv, 'calendar', 'standard'))
                print('First/last time as stored:', str(dates[0]), str(dates[1]))
                for attr in ('history', 'description', 'time_zone', 'timezone'):
                    if attr in ds.ncattrs():
                        print(attr, ds.getncattr(attr))

def sounding_times(year, cfg):
    offset = cfg['local_utc_offset_hours']
    return [t for t in hourly_times(year, cfg) if (t - timedelta(hours=offset)).hour in (0, 6, 12, 18)]

def verify_annual_netcdf(cfg, kind, year, grid_rows):
    """Cheap annual calendar/grid preflight before launching thousands of tasks.

    Step 1 scans all required values; this verifies its completion marker and
    actual coordinates/calendar. Per-grid conversion still checks values.
    """
    import netCDF4 as nc
    path = input_paths(cfg, kind, year)[0]
    manifest_path = path.with_suffix('.manifest.json')
    require(manifest_path.is_file(), f'Missing Step-1 manifest: {manifest_path}')
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get('required_input_complete') is True, f'{path}: Step 1 reports incomplete required input; fix raw gaps before annual AERMET.')
    stat = path.stat()
    require(manifest.get('output_stat') == {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}, f'{path}: file changed since Step-1 completion; rerun/verify preprocessing.')
    with nc.Dataset(path) as ds:
        require(str(getattr(ds, 'aermet_required_input_complete', '')).lower() == 'true', f'{path}: incomplete annual input')
        require(str(getattr(ds, 'time_basis', '')) == 'UTC', f'{path}: expected UTC')
        ids = [str(int(x)) for x in ds.variables['grid_ID'][:]]
        require(len(ids) == len(set(ids)), f'{path}: duplicate grid IDs')
        index = {gid: i for i, gid in enumerate(ids)}
        for row in grid_rows:
            require(row['grid_ID'] in index, f"{path}: missing grid {row['grid_ID']}")
            i = index[row['grid_ID']]
            require(abs(float(ds.variables['latitude'][i]) - float(row['center_lat'])) < 1e-06 and abs(float(ds.variables['longitude'][i]) - float(row['center_lon'])) < 1e-06, f"{path}: coordinate mismatch for {row['grid_ID']}")
        tv = ds.variables['time']
        decoded = nc.num2date(tv[:], tv.units, calendar=getattr(tv, 'calendar', 'standard'))
        stored = [datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in decoded]
        require(all((t.minute == t.second == 0 for t in stored)), f'{path}: off-hour time')
        require(all((b > a for a, b in zip(stored, stored[1:]))), f'{path}: unsorted/duplicate time')
        local = {t + timedelta(hours=local_utc_offset_hours) for t in stored}
        expected = hourly_times(year, cfg) if kind == 'surface' else sounding_times(year, cfg)
        require(set(expected) <= local, f'{path}: annual calendar incomplete')

def completed_matches(folder, cfg):
    done = folder / 'DONE.json'
    if not folder.exists():
        return False
    require(resume_completed and done.is_file(), f'Existing incomplete folder retained: {folder}. Move it aside explicitly or use a new output_dir before retrying.')
    saved = json.loads(done.read_text())
    require(saved['configuration'] == cfg, f'Configuration changed for {folder}; choose a new output_dir.')
    changed_script = saved.get('script_sha256') != script_hash()
    require(not changed_script or recheck_previous_script_outputs, f'Script changed for {folder}; output recheck is disabled.')
    for kind in ('surface', 'upper'):
        for info in saved['source_files'][kind]:
            current = file_info(info['path'])
            require(all((current[k] == info[k] for k in ('path', 'size', 'mtime_ns'))), f'Meteorology changed since completed case: {folder}')
    for info in [saved['aersurface'], saved['executable']] + list(saved.get('output_files', {}).values()):
        require(file_info(info['path']) == info, f"Input/output/executable changed: {info['path']}")
    require(len(saved.get('output_files', {})) == 2, f'Missing saved output provenance: {folder}')
    elev, elev_source = elevation_for(cfg, {'grid_ID': str(saved['grid_id'])}) if not cfg.get('elevation_column') else (saved['elevation_m'], None)
    require(elev == saved['elevation_m'], f'Terrain elevation changed: {folder}')
    if changed_script:
        marker = folder / 'PRESSURE_FIX_RESUME_RECHECK.json'
        if marker.is_file():
            check = json.loads(marker.read_text())
            require(check.get('output_files') == saved['output_files'], f'Existing migration check does not match current script: {folder}')
            if check.get('new_script_sha256') != script_hash():
                audit_outputs(folder, int(saved['year']), cfg['maximum_noncalm_invalid_hours'], cfg)
        else:
            audit = audit_outputs(folder, int(saved['year']), cfg['maximum_noncalm_invalid_hours'], cfg)
            dump(marker, {'old_script_sha256': saved.get('script_sha256'), 'new_script_sha256': script_hash(), 'output_files': saved['output_files'], 'output_audit': audit, 'reason': 'pressure-precision fix; already successful case kept unchanged'})
    return True

def archive_known_pressure_failure(folder, cfg):
    if not folder.exists() or (folder / 'DONE.json').exists():
        return
    manifest = folder / 'manifest.json'
    if not manifest.is_file():
        if archive_orphan_folders_from_stopped_job:
            root = Path(cfg['output_root']).resolve()
            resolved = folder.resolve()
            require(not folder.is_symlink() and resolved.parent.parent == root and re.fullmatch('20\\d\\d', resolved.parent.name) and re.fullmatch('Grid_\\d+', resolved.name), f'Unexpected orphan folder: {folder}')
            archive = folder.with_name(folder.name + '_orphan_stopped_job_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
            require(archive.resolve().parent == resolved.parent and (not archive.exists()), 'Conflicting orphan archive')
            folder.rename(archive)
            logging.info('Preserved manifest-less folder in %s', archive)
        return
    saved = json.loads(manifest.read_text())
    pressure_failure = archive_pressure_rounding_failures and saved.get('state') == 'failed' and (saved.get('error') == 'Nondecreasing rounded profile pressure')
    unfinished = archive_unfinished_cases_from_stopped_job and saved.get('state') in {'started', 'prepared_only'}
    messages = folder / 'stage2.mes'
    stage2_elevation_failure = False
    if archive_confirmed_upperair_elevation_failures and saved.get('state') == 'failed' and (saved.get('error') == 'stage2: no unequivocal successful completion; see report/logs') and messages.is_file():
        text = messages.read_text(errors='replace').upper()
        stage2_elevation_failure = any(('UPPERAIR' in line and 'E05' in line and ('INVALID FORMAT FOR STATION ELEVATION:' in line) and ('*******' in line) for line in text.splitlines()))
    stage1_missing_elevation = False
    messages1 = folder / 'AERMET_MSG.MES'
    if archive_confirmed_upperair_elevation_failures and saved.get('state') == 'failed' and (saved.get('error') == 'stage1: no unequivocal successful completion; see report/logs') and messages1.is_file():
        stage1_missing_elevation = any(('UPPERAIR' in line and 'E05' in line and ('STATION ELEVATION NOT ENTERED' in line) for line in messages1.read_text(errors='replace').upper().splitlines()))
    if not (pressure_failure or unfinished or stage2_elevation_failure or stage1_missing_elevation):
        return
    require(saved.get('configuration') == cfg, f'Failed-case configuration changed: {folder}')
    root = Path(cfg['output_root']).resolve()
    resolved = folder.resolve()
    require(resolved.parent.parent == root and resolved.name == f"Grid_{saved['grid_id']}" and (resolved.parent.name == str(saved['year'])), f'Unexpected archive path: {resolved}')
    tag = '_failed_pressure_' if pressure_failure else '_failed_UA_elevation_' if stage2_elevation_failure or stage1_missing_elevation else '_unfinished_stopped_job_'
    archive = folder.with_name(folder.name + tag + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    require(archive.resolve().parent == resolved.parent and (not archive.exists()), 'Unsafe/conflicting archive target')
    folder.rename(archive)
    logging.info('Preserved failed case in %s; retrying original case path', archive)

def production_task(cfg, row, year):
    """Always return a status so one failure does not hide other case results."""
    folder = Path(cfg['output_root']) / str(year) / f"Grid_{row['grid_ID']}"
    try:
        archive_known_pressure_failure(folder, cfg)
        if completed_matches(folder, cfg):
            return {'grid_id': row['grid_ID'], 'year': year, 'state': 'skipped_verified_complete', 'folder': str(folder)}
        return process_grid(cfg, row, year, prepare_only)
    except Exception as exc:
        logging.error('FAILED Grid %s Year %s: %s', row['grid_ID'], year, exc)
        return {'grid_id': row['grid_ID'], 'year': year, 'state': 'failed', 'folder': str(folder), 'error': str(exc)}

def append_summary(path, result):
    with path.open('a', encoding='utf-8') as f:
        f.write(json.dumps(result, default=str) + '\n')
    logging.info('Grid %s Year %s: %s', result['grid_id'], result['year'], result['state'])

def main():
    cfg = build_configuration()
    require(period_start is None and period_end is None, 'This is the annual runner; both period limits must be None.')
    require(surface_time_basis == upper_time_basis == 'UTC' and local_utc_offset_hours == -5, 'This validated rebuild/run pair requires UTC inputs, fixed EST offset -5, no DST.')
    require(max_workers >= 1 and years and (len(set(years)) == len(years)), 'Check workers/years')
    allocated = int(__import__('os').environ.get('SLURM_CPUS_PER_TASK', str(max_workers)))
    require(max_workers <= allocated, 'max_workers exceeds CPUs allocated by the job.')
    with Path(grid_csv_file).open(newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r['grid_ID'] = str(int(float(r['grid_ID'])))
    require(len({r['grid_ID'] for r in rows}) == len(rows), 'Duplicate grid IDs')
    if grid_ids is not None:
        wanted = {str(int(g)) for g in grid_ids}
        rows = [r for r in rows if r['grid_ID'] in wanted]
        require({r['grid_ID'] for r in rows} == wanted, 'Some selected grids not found')
    require(rows, 'No grids selected')
    if inspect_only:
        for year in years:
            inspect_netcdf_inputs(cfg, year)
        return 0
    if not prepare_only:
        require(Path(aermet_executable).is_file(), 'AERMET executable not found')
    errors = []
    for year in years:
        for kind in ('surface', 'upper'):
            try:
                verify_annual_netcdf(cfg, kind, year, rows)
            except Exception as exc:
                errors.append(str(exc))
        for r in rows:
            try:
                require(Path(cfg['aersurface_pattern'].format(year=year, grid_id=r['grid_ID'])).is_file(), f"Missing AERSURFACE: year={year} grid={r['grid_ID']}")
                elevation_for(cfg, r)
            except Exception as exc:
                errors.append(str(exc))
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    dump(output / f'PREFLIGHT_{run_id}.json', {'configuration': cfg, 'grids': len(rows), 'years': list(years), 'errors': errors})
    require(not errors, f'Annual preflight found {len(errors)} problems; see PREFLIGHT_{run_id}.json. First: {errors[:3]}')
    summary = output / f'RUN_SUMMARY_{run_id}.jsonl'
    failures = 0
    logging.info('Annual AERMET: %s grids x %s years; workers=%s; output=%s', len(rows), len(years), max_workers, output)
    for year in years:
        if max_workers == 1:
            for r in rows:
                result = production_task(cfg, r, year)
                append_summary(summary, result)
                failures += result['state'] == 'failed'
        else:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                futures = {executor.submit(production_task, cfg, r, year): r['grid_ID'] for r in rows}
                for future in as_completed(futures):
                    try:
                        result = future.result()
                    except Exception as exc:
                        result = {'grid_id': futures[future], 'year': year, 'state': 'failed', 'error': str(exc)}
                    append_summary(summary, result)
                    failures += result['state'] == 'failed'
    logging.info('Finished annual run; failed=%s; summary=%s', failures, summary)
    return 1 if failures else 0
if __name__ == '__main__':
    result_code = main()
    if result_code:
        raise RuntimeError('Some AERMET cases failed. Read RUN_SUMMARY and per-grid failure.txt. Existing cases were retained.')
