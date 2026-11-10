"""Complete documented surface-field gaps using the study substitution rules."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
from pathlib import Path
from datetime import datetime, timedelta
import csv
import hashlib
import json
import logging
import math
import numpy as np
import netCDF4 as nc
input_folder = Path(_configured_path('complete_surface_fields.input_folder_1'))
output_folder = Path(_configured_path('complete_surface_fields.output_folder_2'))
diagnostic_file = Path(_configured_path('surface_inspection_report'))
years = _setting('complete_surface_fields.years')
maximum_radiation_donor_km = 30.0
maximum_negative_tp_m = 5e-06
fill_boundary_hours = True
resume_completed = True
FIELDS = ('t2m', 'd2m', 'u10', 'v10', 'sp', 'tp', 'tcc', 'ssrd', 'strd')
REQUIRED = FIELDS[:-1]
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

def require(ok, message):
    if not ok:
        raise ValueError(message)

def dump(path, data):
    path.write_text(json.dumps(data, indent=2, default=str) + '\n', encoding='utf-8')

def stat_info(path):
    s = path.stat()
    return {'size': s.st_size, 'mtime_ns': s.st_mtime_ns}

def data_array(variable, key):
    return np.ma.asarray(variable[key]).astype(float).filled(np.nan)

def nearest_donor(index, good, lat, lon):
    a = np.radians(lat[good] - lat[index])
    b = np.radians((lon[good] - lon[index] + 180) % 360 - 180)
    h = np.sin(a / 2) ** 2 + np.cos(np.radians(lat[index])) * np.cos(np.radians(lat[good])) * np.sin(b / 2) ** 2
    distances = 12742 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))
    j = int(np.argmin(distances))
    return (int(good[j]), float(distances[j]))

def completion_year(year, report, diagnostic_sha):
    name = f'{year}_surface_variables_yearly.nc'
    source = input_folder / name
    destination = output_folder / name
    manifest_path = destination.with_suffix('.manifest.json')
    original_manifest = source.with_suffix('.manifest.json')
    settings = {'radiation_max_km': maximum_radiation_donor_km, 'tp_negative_cap_m': maximum_negative_tp_m, 'fill_boundary_hours': fill_boundary_hours, 'diagnostic_sha256': diagnostic_sha, 'algorithm': 'same_hour_original_nearest_ssrd_verified_small_tp_zero_boundary03UTC_v1'}
    if destination.exists() or manifest_path.exists():
        require(resume_completed and destination.is_file() and manifest_path.is_file(), f'Existing output retained: {destination}')
        saved = json.loads(manifest_path.read_text())
        require(saved.get('required_input_complete') is True and saved.get('completion_settings') == settings and (saved.get('completion_input_stat') == stat_info(source)) and (saved.get('output_stat') == stat_info(destination)) and (saved.get('completion_original_manifest_sha256') == hashlib.sha256(original_manifest.read_bytes()).hexdigest()), f'Existing completion changed/incomplete: {destination}; use a new output folder.')
        return {'year': year, 'state': 'skipped_verified_complete', 'output': str(destination)}
    partial = destination.with_name(name + '.partial')
    audit_file = destination.with_suffix('.completion.csv')
    require(not partial.exists() and (not audit_file.exists()), f'Existing partial/audit retained: {partial}')
    manifest = json.loads(original_manifest.read_text())
    events = {}
    for event in report['events']:
        require('diagnostic_error' not in event, f'Diagnostic error: {event}')
        key = (event['variable'], event['time_UTC'])
        require(key not in events, f'Duplicate diagnostic event: {key}')
        require(len(event['sample_grid_IDs']) == event['affected_grids'], f'Diagnostic sample truncated at {key}; rerun diagnosis with larger sample limit.')
        events[key] = event
    changed = 0
    counts = {'tp_zero': 0, 'ssrd_spatial_nearest': 0, 'boundary_nearest_time': 0}
    with nc.Dataset(source) as src:
        tv = src.variables['time']
        calendar = getattr(tv, 'calendar', 'standard')
        decoded = nc.num2date(tv[:], tv.units, calendar=calendar)
        old_times = [datetime(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in decoded]
        require(all((b > a for a, b in zip(old_times, old_times[1:]))), f'Unsorted/duplicate times in {source}')
        old_index = {t: i for i, t in enumerate(old_times)}
        required = [datetime(year, 1, 1) + timedelta(hours=h) for h in range(6, int((datetime(year + 1, 1, 1) - datetime(year, 1, 1)).total_seconds() / 3600) + 6)]
        required_set = set(required)
        missing = required_set - set(old_times)
        authorized = {datetime(2024, 1, 1, 4), datetime(2024, 1, 1, 5)} if year == 2023 and fill_boundary_hours else set()
        require(missing <= authorized, f'Unapproved missing timestamps in {year}: {sorted(missing)}')
        require(not missing or datetime(2024, 1, 1, 3) in old_index, 'Missing boundary donor 03 UTC')
        times = sorted(set(old_times) | missing)
        index = {t: i for i, t in enumerate(times)}
        mapping = np.array([old_index[t] if t in old_index else old_index[datetime(2024, 1, 1, 3)] for t in times])
        lat = data_array(src.variables['latitude'], slice(None))
        lon = data_array(src.variables['longitude'], slice(None))
        ids = [str(int(i)) for i in src.variables['grid_ID'][:]]
        require(len(set(ids)) == len(ids), 'Duplicate grid IDs')
        fields = ['year', 'variable', 'time_UTC', 'grid_ID', 'action', 'old_value', 'new_value', 'donor_grid_ID', 'donor_time_UTC', 'distance_km']
        with audit_file.open('x', newline='', encoding='utf-8') as audit, nc.Dataset(partial, 'w', format='NETCDF4') as out:
            writer = csv.DictWriter(audit, fieldnames=fields)
            writer.writeheader()

            def record(variable, t, j, action, old, new, donor=None, donor_time=None, km=None):
                nonlocal changed
                changed += 1
                counts[action] += 1
                writer.writerow(dict(year=year, variable=variable, time_UTC=str(t), grid_ID=ids[j], action=action, old_value=old, new_value=float(new), donor_grid_ID=ids[donor] if donor is not None else '', donor_time_UTC=str(donor_time) if donor_time is not None else '', distance_km=km if km is not None else ''))
            out.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
            out.product_status = 'completion_in_progress'
            for dim, obj in src.dimensions.items():
                out.createDimension(dim, len(times) if dim == 'time' else len(obj))
            for vn, var in src.variables.items():
                attrs = {k: var.getncattr(k) for k in var.ncattrs() if k != '_FillValue'}
                kwargs = {'zlib': True, 'complevel': 2} if 'time' in var.dimensions else {}
                if '_FillValue' in var.ncattrs():
                    kwargs['fill_value'] = var.getncattr('_FillValue')
                dest = out.createVariable(vn, var.datatype, var.dimensions, **kwargs)
                dest.setncatts(attrs)
                if vn == 'time':
                    dest[:] = nc.date2num(times, tv.units, calendar=calendar)
                elif vn == 'time_bnds':
                    values = nc.date2num(times, tv.units, calendar=calendar)
                    require(tv.units.startswith('hours since'), 'Expected hourly coordinate units')
                    dest[:] = np.column_stack([values - 1, values])
                elif 'time' not in var.dimensions:
                    dest[:] = var[:]
                else:
                    require(var.dimensions == ('time', 'location'), f'Unsupported timed variable: {vn} {var.dimensions}')
                    for lo in range(0, len(times), 168):
                        hi = min(lo + 168, len(times))
                        dest[lo:hi, :] = var[mapping[lo:hi], :]
            visited = set()
            for name in REQUIRED:
                require(name in out.variables, f'Missing variable {name}')
                v = out.variables[name]
                for lo in range(0, len(times), 168):
                    hi = min(lo + 168, len(times))
                    block = data_array(v, (slice(lo, hi), slice(None)))
                    for row, j in np.argwhere(~np.isfinite(block)):
                        i = lo + int(row)
                        j = int(j)
                        t = times[i]
                        if t not in required_set:
                            continue
                        require(t not in missing, f'Boundary donor itself has invalid {name}: {ids[j]}')
                        key = (name, str(t))
                        event = events.get(key)
                        require(event is not None, f'Undiagnosed missing {name} {t} grid {ids[j]}')
                        require(ids[j] in event['sample_grid_IDs'], f'Grid not in diagnostic: {key} {ids[j]}')
                        visited.add((name, str(t), ids[j]))
                        if name == 'tp':
                            k = event['sample_grid_IDs'].index(ids[j])
                            now = event['current_interpolated_values'][k]
                            previous = event['previous_interpolated_values'][k]
                            require(now is not None and previous is not None, 'Raw precipitation is missing; cannot classify as negative-roundoff correction')
                            delta = now - previous
                            require(-maximum_negative_tp_m <= delta < 0, f'Negative precipitation exceeds correction cap: {delta}')
                            block[row, j] = 0
                            record(name, t, j, 'tp_zero', delta, 0)
                        elif name == 'ssrd':
                            original = data_array(src.variables[name], (old_index[t], slice(None)))
                            good = np.flatnonzero(np.isfinite(original) & (original >= 0))
                            require(len(good) > 0, f'No original same-hour radiation donor: {t}')
                            donor, km = nearest_donor(j, good, lat, lon)
                            require(km <= maximum_radiation_donor_km, f'Nearest radiation donor too far: {km:.2f} km at {t}')
                            value = float(original[donor])
                            block[row, j] = value
                            record(name, t, j, 'ssrd_spatial_nearest', 'missing', value, donor, t, km)
                        else:
                            raise ValueError(f'No approved completion for missing {name}')
                    v[lo:hi, :] = block
            expected = {(name, t, gid) for (name, t), event in events.items() for gid in event['sample_grid_IDs']}
            require(visited == expected, f'Diagnostic does not match existing NetCDF gaps: found={len(visited)}, expected={len(expected)}')
            for t in sorted(missing):
                donor_time = datetime(2024, 1, 1, 3)
                for name in FIELDS:
                    vals = data_array(out.variables[name], (index[t], slice(None)))
                    for j, value in enumerate(vals):
                        require(np.isfinite(value), f'Boundary donor not finite: {name} grid {ids[j]}')
                        record(name, t, j, 'boundary_nearest_time', 'missing', value, j, donor_time, 0)
            for name in REQUIRED:
                for lo in range(0, len(times), 168):
                    values = data_array(out.variables[name], (slice(lo, lo + 168), slice(None)))
                    selected = [t in required_set for t in times[lo:lo + 168]]
                    require(np.isfinite(values[selected, :]).all(), f'Remaining required missing {name}')
            out.product_status = 'numerically_complete_with_documented_estimates'
            out.aermet_required_input_complete = 'true'
            out.missing_required_UTC_timestamps = '[]'
            out.completion_audit_file = audit_file.name
            out.completion_estimated_boundary_hours = json.dumps([str(t) for t in sorted(missing)])
            out.completion_method = json.dumps(settings)
            out.history = str(getattr(src, 'history', '')) + '; diagnosed gaps completioned in a new file; see completion CSV'
    partial.rename(destination)
    manifest['pre_completion_diagnostics'] = {k: manifest.get(k) for k in ['variable_diagnostics', 'missing_required_UTC_times', 'missing_required_raw_times_by_variable']}
    manifest.update(state='completed_with_documented_completions', output=str(destination), required_input_complete=True, missing_required_UTC_times=[], completion_settings=settings, completion_counts=counts, completion_audit_file=str(audit_file), completion_input_stat=stat_info(source), completion_original_manifest_sha256=hashlib.sha256(original_manifest.read_bytes()).hexdigest(), output_stat=stat_info(destination), readiness_note='Required values now finite, including documented estimates. Raw retrieval gaps have not been eliminated.')
    dump(manifest_path, manifest)
    logging.info('COMPLETIONED %s: %s', year, counts)
    return {'year': year, 'state': 'completed_with_documented_completions', 'output': str(destination), 'counts': counts}

def main():
    require(input_folder.resolve() != output_folder.resolve(), 'Input/output folders must differ')
    report = json.loads(diagnostic_file.read_text())
    lookup = {int(r['year']): r for r in report}
    sha = hashlib.sha256(diagnostic_file.read_bytes()).hexdigest()
    output_folder.mkdir(parents=True, exist_ok=True)
    results = []
    for year in years:
        try:
            require(year in lookup, f'No diagnostic report for {year}')
            results.append(completion_year(year, lookup[year], sha))
        except Exception as exc:
            logging.exception('FAILED year %s; originals retained', year)
            results.append({'year': year, 'state': 'failed', 'error': str(exc)})
        dump(output_folder / 'SURFACE_COMPLETION_SUMMARY.json', results)
    require(not any((r['state'] == 'failed' for r in results)), 'Some completions failed; read SURFACE_COMPLETION_SUMMARY.json. Existing/partial outputs were retained.')
    print('Surface completions complete. Share SURFACE_COMPLETION_SUMMARY.json before launching full AERMET.')
if __name__ == '__main__':
    main()
