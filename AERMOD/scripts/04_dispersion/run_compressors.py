"""Run AERMOD for compressor emissions."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
import os, shutil, math, re, logging, subprocess
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import utm
AERMOD_EXECUTABLE = _configured_path('run_compressors.aermod_executable_1')
METEOROLOGICAL_DATA_DIR_TEMPLATE = _configured_path('run_compressors.meteorological_data_dir_template_2')
GRID_CSV_FILE = _configured_path('run_compressors.grid_csv_file_3')
AERMAP_GRIDWISE_DIR = _configured_path('run_compressors.aermap_gridwise_dir_4')
YEARS = _setting('run_compressors.years')
POLLUTANTS = _setting('run_compressors.pollutants')
SOURCE_IDS = _setting('run_compressors.source_ids')
MAX_WORKERS = _setting('run_compressors.max_workers')
EXPECTED_AERMOD_VERSION = '24142'
MODEL_TIMEOUT_SECONDS = 21600
LINK_MET_FILES = True
ALLOW_DOCUMENTED_NEUTRAL_EXCEPTIONS = True
MASS_GRAMS_PER_TON = 1000000.0
WORKFLOW = 'compressors'
OUTPUT_BASE_DIR_TEMPLATE = _configured_path('run_compressors.output_base_dir_template_5')
LOG_DIRECTORY = _configured_path('run_compressors.log_directory_6')
INVENTORY_TEMPLATE = _configured_path('run_compressors.inventory_template_7')

def sanitize_identifier(identifier):
    """Remove unwanted characters from an identifier."""
    return re.sub('\\W+', '', str(identifier))

def standardize_facility_id(x):
    """
    Convert the facility ID to a standardized string.
    If the value is numeric (or can be converted to a float) and represents an integer,
    return the integer value as a string.
    Otherwise, just return the trimmed string.
    """
    s = str(x).strip()
    try:
        num = float(s)
        if num.is_integer():
            return str(int(num))
        else:
            return str(num)
    except ValueError:
        return s

def _write_control_file(facility_id, year, stack_params, hourly_emis_file, output_dir, surffile_path, proffile_path, combined_receptor_file, elevation, pollutant):
    """
    Generate the AERMOD .inp file with all required sections.
    """
    inp_file = os.path.join(output_dir, f'aermod_{facility_id}.inp')
    error_file = os.path.join(output_dir, f'error_{facility_id}.err')
    summfile_path = os.path.join(output_dir, f'aermod_{facility_id}.sum')
    source_id = f'{facility_id}'[:8]
    try:
        with open(inp_file, 'w') as f:
            f.write('CO STARTING\n')
            f.write(f'   TITLEONE  AERMOD Run for {facility_id}\n')
            if pollutant.upper() == 'NO2':
                f.write('   MODELOPT  CONC ALPHA WETDPLT TTRM2 ARM2\n')
                f.write('   OZONEVAL  40.0\n')
                f.write('   NO2STACK  0.10\n')
            elif pollutant.upper() == 'PM2.5':
                f.write('   MODELOPT  CONC ALPHA WETDPLT DRYDPLT\n')
            elif pollutant.upper() == 'SO2':
                f.write('   MODELOPT  CONC ALPHA WETDPLT\n')
            else:
                f.write('   MODELOPT  CONC\n')
            f.write('   AVERTIME  24\n')
            f.write(f'   POLLUTID  {pollutant}\n')
            f.write('   RUNORNOT  RUN\n')
            f.write(f'   ERRORFIL  "{error_file}"\n')
            f.write('CO FINISHED\n\n')
            f.write('SO STARTING\n')
            f.write('   ELEVUNIT  METERS\n')
            if pollutant.upper() == 'NO2':
                f.write(f'   GASDEPOS  {source_id} 0.1361 1.4E-5 12000 8444\n')
            elif pollutant.upper() == 'PM2.5':
                f.write(f'   PARTDIAM {source_id} 0.1\n')
                f.write(f'   PARTDENS {source_id} 1.0\n')
                f.write(f'   MASSFRAX {source_id} 1.0\n')
            elif pollutant.upper() == 'SO2':
                f.write(f'   GASDEPOS  {source_id} 0.1122 1.83E-5 732 72\n')
            f.write('   SRCGROUP  ALL\n')
            f.write('SO FINISHED\n\n')
            if combined_receptor_file and os.path.exists(combined_receptor_file):
                f.write('RE STARTING\n')
                f.write(f'   INCLUDED  "{combined_receptor_file}"\n')
                f.write('RE FINISHED\n\n')
            f.write('ME STARTING\n')
            f.write(f'   SURFFILE "{surffile_path}"\n')
            f.write(f'   PROFFILE "{proffile_path}"\n')
            f.write(f'   SURFDATA 99999 {year} STUDY_AREA\n')
            f.write(f'   UAIRDATA 99999 {year} STUDY_AREA\n')
            f.write('   PROFBASE ' + str(json.loads((Path(output_dir) / 'met_context.json').read_text())['elevation_m']) + '\n')
            f.write(f'   STARTEND {year} 01 01 {year} 12 31\n')
            f.write('ME FINISHED\n\n')
            f.write('OU STARTING\n')
            f.write(f'   SUMMFILE "{summfile_path}"\n')
            f.write('   DAYTABLE 24\n')
            f.write('OU FINISHED\n\n')
        logging.info(f'AERMOD input file created: {inp_file}')
        return inp_file
    except Exception as e:
        logging.error(f'Failed to create AERMOD input file: {e}')
        return None
from pathlib import Path
from functools import lru_cache
from decimal import Decimal, InvalidOperation
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import calendar, csv, hashlib, json, traceback

def require(ok, message):
    if not ok:
        raise ValueError(message)

def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=str, allow_nan=False) + '\n', encoding='utf-8')

def file_record(path):
    p = Path(path)
    s = p.stat()
    return {'path': str(p.resolve()), 'size': s.st_size, 'mtime_ns': s.st_mtime_ns}

def log_and_print(message, level='info', silent=False):
    if level in ('error', 'critical'):
        print(message, flush=True)

def ensure_directory_exists(path):
    Path(path).mkdir(parents=True, exist_ok=True)

def identifier(raw):
    value = str(raw).strip()
    require(value and value.lower() not in {'nan', 'none', ''}, 'Missing source ID')
    try:
        n = Decimal(value)
        if n.is_finite() and n == n.to_integral_value():
            value = str(int(n))
    except InvalidOperation:
        pass
    if WORKFLOW == 'activity':
        value = value.zfill(10)
    require(re.fullmatch('[A-Za-z0-9_.-]+', value) is not None, f'Unsafe source ID: {value}')
    return value

def normalize_columns(df):
    df.columns = [str(c).strip().lower() for c in df.columns]
    return df

@lru_cache(maxsize=2)
def table(path):
    return normalize_columns(pd.read_csv(path, dtype=str))

def source_row(path, source):
    df = table(str(path))
    key = id_column()
    require(key in df.columns, f'{path}: missing {key}')
    selected = df[df[key].map(identifier) == source]
    require(len(selected) == 1, f'{path}: expected one row for {source}; found {len(selected)}')
    return selected.iloc[0]

def numeric(value, default=None):
    if pd.isna(value) or str(value).strip().lower() in {'', 'nan', 'none'}:
        require(default is not None, 'Required numeric value is missing')
        return default
    result = float(value)
    require(math.isfinite(result), 'Nonfinite numeric input')
    return result

@lru_cache(maxsize=1)
def grids():
    df = table(GRID_CSV_FILE).copy()
    df['grid_id'] = df['grid_id'].map(lambda v: str(int(Decimal(str(v)))))
    require(not df['grid_id'].duplicated().any(), 'Duplicate grid IDs')
    for c in ['center_lat', 'center_lon']:
        df[c] = pd.to_numeric(df[c], errors='raise')
    return df

def distances(lat, lon):
    df = grids()
    la = df['center_lat'].to_numpy(float)
    lo = df['center_lon'].to_numpy(float)
    a = np.sin(np.radians(la - lat) / 2) ** 2 + math.cos(math.radians(lat)) * np.cos(np.radians(la)) * np.sin(np.radians(lo - lon) / 2) ** 2
    return 12742 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

def _prepare_meteorology(folder, gid, year):
    source = Path(METEOROLOGICAL_DATA_DIR_TEMPLATE.format(year=year, grid_ID=gid))
    manifest = json.loads((source / 'manifest.json').read_text())
    audit = json.loads((source / 'output_audit.json').read_text())
    require(str(int(manifest['grid_id'])) == gid and int(manifest['year']) == year, 'Meteorology manifest grid/year mismatch')
    require(audit['hours'] == (8784 if calendar.isleap(year) else 8760), 'Incorrect annual meteorology count')
    bad = int(audit['noncalm_invalid_hours'])
    exceptions = []
    if bad:
        require(ALLOW_DOCUMENTED_NEUTRAL_EXCEPTIONS and bad == 1 and (len(audit['invalid_times']) == 1), f'Unapproved missing meteorology: Grid {gid} year {year}')
        wanted = set(audit['invalid_times'])
        for line in (source / 'output_stage2.SFC').read_text().splitlines()[1:]:
            p = line.split()
            if not p:
                continue
            y = int(p[0])
            y = y + 2000 if y < 100 else y
            stamp = datetime(y, int(p[1]), int(p[2])) + timedelta(hours=int(p[4]))
            if str(stamp) in wanted:
                require(float(p[5]) == 0 and float(p[7]) == -9 and (float(p[9]) == -999) and (float(p[11]) == -8888), 'Missing hour differs from reviewed near-neutral pattern')
                exceptions.append({'time_EST': str(stamp), 'SFC_record': line})
        require(len(exceptions) == bad, 'Exception record not found')
    else:
        require((source / 'DONE.json').is_file() and manifest.get('state') == 'completed_pending_scientific_review', 'Meteorology has no completed status')
    for report in ['AERMET_REPORT.RPT', 'stage2.rpt']:
        text = (source / report).read_text(errors='replace').upper()
        require('AERMET FINISHED SUCCESSFULLY' in text and 'UN-SUCCESSFULLY' not in text, 'AERMET stage not successful')
    for name in ['output_stage2.SFC', 'output_stage2.PFL']:
        require((source / name).stat().st_size > 0, 'Empty meteorological file')
        if LINK_MET_FILES:
            (folder / name).symlink_to((source / name).resolve())
        else:
            shutil.copy2(source / name, folder / name)
    info = {'grid_id': gid, 'year': year, 'current_folder': str(source), 'elevation_m': numeric(manifest['elevation_m']), 'exceptions': exceptions, 'files': {n: file_record(source / n) for n in ['output_stage2.SFC', 'output_stage2.PFL']}}
    dump(folder / 'met_context.json', info)
    return info

def model_run(inp, folder, year):
    inp = Path(inp)
    text = inp.read_text()
    context = json.loads((folder / 'met_context.json').read_text())
    text = re.sub('(?m)^\\s*PROFBASE\\s+.*$', f"   PROFBASE {context['elevation_m']:.3f}", text)
    text = re.sub('(?m)^\\s*ERRORFIL\\s+.*\\n', '', text)
    text = re.sub('(?m)^\\s*SUMMFILE\\s+.*\\n', '', text)
    text = text.replace('CO FINISHED', '   ERRORFIL "model.err"\nCO FINISHED')
    text = text.replace('OU FINISHED', '   SUMMFILE "model.sum"\nOU FINISHED')
    inp.write_text(text, encoding='ascii')
    output = inp.with_suffix('.out')
    with (folder / 'aermod_run.log').open('wb') as log:
        result = subprocess.run([AERMOD_EXECUTABLE, inp.name, output.name], cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=MODEL_TIMEOUT_SECONDS, check=False)
    version = None
    fatal = None
    missing = None
    units = None
    with output.open(errors='replace') as f:
        for line in f:
            m = re.search('AERMOD\\s*-\\s*VERSION\\s+(\\d+)', line, re.I)
            if m:
                version = m[1]
            m = re.search('A Total of\\s+(\\d+)\\s+Fatal Error', line, re.I)
            if m:
                fatal = int(m[1])
            m = re.search('A Total of\\s+(\\d+)\\s+Missing Hours', line, re.I)
            if m:
                missing = int(m[1])
            if 'Output Units' in line:
                units = line.strip()
    require(result.returncode == 0 and version == EXPECTED_AERMOD_VERSION and (fatal == 0), f'AERMOD not successful: return={result.returncode}; version={version}; fatal={fatal}')
    require(missing is not None and missing == len(context['exceptions']), f"AERMOD missing-hour count {missing} differs from reviewed meteorological exceptions {len(context['exceptions'])}")
    return (output, {'version': version, 'fatal_errors': fatal, 'missing_met_hours_reported': missing, 'units': units})

def daily_csv(output, year, receptors):
    n = 366 if calendar.isleap(year) else 365
    values = {xy: np.full(n, np.nan) for xy in receptors}
    current = None
    group = None
    recording = False
    with output.open(errors='replace') as f:
        for line in f:
            if 'CONCURRENT' in line:
                current = None
                group = None
                recording = False
                m = re.search('CONCURRENT\\s+24\\s*(?:-?HR)?\\s+(?:AVERAGE\\s+)?CONC(?:ENTRATION)?\\s+VALUES.*HOUR\\s+(\\d+)\\s+FOR DAY\\s+(\\d+)\\s+OF\\s+(\\d{4})', line)
                if m and int(m[1]) == 24 and (int(m[3]) == year):
                    current = int(m[2]) - 1
                    require(0 <= current < n, 'Invalid Julian day')
            if 'FOR SOURCE GROUP:' in line:
                group = line.split('FOR SOURCE GROUP:', 1)[1].strip().split()[0]
            if '*** DISCRETE CARTESIAN RECEPTOR POINTS ***' in line:
                recording = current is not None and group == 'ALL'
                continue
            if '***' in line or '\x0c' in line:
                recording = False
            if not recording:
                continue
            p = line.split()
            if len(p) not in (3, 6):
                continue
            try:
                numbers = [float(v.replace('D', 'E')) for v in p]
            except ValueError:
                continue
            for i in range(0, len(numbers), 3):
                x, y, c = numbers[i:i + 3]
                xy = (round(x, 2), round(y, 2))
                require(xy in values and math.isfinite(c) and (c >= 0), 'Unexpected receptor/value in daily output')
                old = values[xy][current]
                require(not np.isfinite(old) or old == c, 'Conflicting repeated receptor/day result')
                values[xy][current] = c
    missing = sum((int(np.isnan(v).sum()) for v in values.values()))
    require(missing == 0, f'{missing} receptor/day results absent; refusing to zero-fill CSV')
    destination = output.with_suffix('.csv')
    with destination.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['utm_easting', 'utm_northing'] + [f'day{d}_conc' for d in range(1, n + 1)])
        for (x, y), v in sorted(values.items()):
            w.writerow([x, y, *v])
    return destination

def validate_saved_daily_csv(path, year, expected_receptors):
    """Read back the saved CSV before allowing deletion of the raw model output."""
    days = 366 if calendar.isleap(year) else 365
    seen = set()
    with Path(path).open(newline='') as handle:
        reader = csv.reader(handle)
        require(next(reader, None) == ['utm_easting', 'utm_northing'] + [f'day{d}_conc' for d in range(1, days + 1)], 'Saved CSV header invalid')
        for row in reader:
            require(len(row) == days + 2, 'Saved CSV row length invalid')
            values = np.asarray(row, dtype=float)
            require(np.isfinite(values).all() and (values[2:] >= 0).all(), 'Saved CSV contains invalid values')
            xy = tuple(values[:2])
            require(xy not in seen, 'Saved CSV has duplicate receptors')
            seen.add(xy)
    require(len(seen) == expected_receptors and len(seen) > 0, 'Saved CSV receptor count invalid')

def compact_success_diagnostics(folder):
    """Compress successful-run diagnostics; preserve their contents for review."""
    import gzip
    import shutil
    for name in ('model.err', 'model.sum', 'aermod_run.log'):
        source = folder / name
        if not source.is_file() or source.is_symlink():
            continue
        target = folder / (name + '.gz')
        temporary = folder / (name + '.gz.tmp')
        with source.open('rb') as inp, gzip.open(temporary, 'wb') as out:
            shutil.copyfileobj(inp, out)
        temporary.replace(target)
        source.unlink()
    (folder / 'AERMOD_STATUS.json').unlink(missing_ok=True)

def consolidated_done(folder):
    """Read a completion record without recreating per-case audit files."""
    import sqlite3
    database = folder.parent / 'audit.sqlite3'
    if not database.is_file():
        return None
    connection = sqlite3.connect(database.as_uri() + '?mode=ro', uri=True, timeout=120)
    try:
        row = connection.execute('SELECT record FROM cases WHERE source=?', (folder.name,)).fetchone()
        return json.loads(row[0]) if row else None
    finally:
        connection.close()

def consolidate_case_audit(folder):
    """Parent process only: commit and verify audits, then remove loose copies.

    One SQLite file per year/pollutant folder stores JSON records and compressed
    diagnostic bytes. SQLite transactions serialize separate job writers.
    Failed cases and model inputs/concentration CSVs are left untouched.
    """
    import sqlite3
    import zlib
    done = folder / 'AERMOD_DONE.json'
    if not done.is_file():
        return
    record = done.read_text()
    require(json.loads(record).get('state') == 'completed', 'Cannot consolidate incomplete case')
    names = ['AERMOD_DONE.json', 'AERMOD_STATUS.json', 'met_context.json', 'model.err', 'model.sum', 'aermod_run.log', 'model.err.gz', 'model.sum.gz', 'aermod_run.log.gz']
    files = {}
    for name in names:
        path = folder / name
        if path.is_file() and (not path.is_symlink()):
            files[name] = path.read_bytes()
    database = folder.parent / 'audit.sqlite3'
    connection = sqlite3.connect(database, timeout=120)
    try:
        connection.execute('PRAGMA synchronous=FULL')
        connection.execute('CREATE TABLE IF NOT EXISTS cases (source TEXT PRIMARY KEY, record TEXT NOT NULL)')
        connection.execute('CREATE TABLE IF NOT EXISTS diagnostics (source TEXT, name TEXT, content BLOB NOT NULL, PRIMARY KEY(source,name))')
        with connection:
            connection.execute('INSERT OR REPLACE INTO cases VALUES (?,?)', (folder.name, record))
            for name, content in files.items():
                connection.execute('INSERT OR REPLACE INTO diagnostics VALUES (?,?,?)', (folder.name, name, zlib.compress(content)))
        require(connection.execute('SELECT record FROM cases WHERE source=?', (folder.name,)).fetchone()[0] == record, 'Audit record verification failed')
        for name, content in files.items():
            saved = connection.execute('SELECT content FROM diagnostics WHERE source=? AND name=?', (folder.name, name)).fetchone()[0]
            require(zlib.decompress(saved) == content, 'Archived diagnostic verification failed')
        for name, content in files.items():
            path = folder / name
            require(not path.is_symlink() and path.read_bytes() == content, 'Audit changed during consolidation; loose copy retained')
            path.unlink()
    finally:
        connection.close()

def _run_model_case(task):
    year, pollutant, source = task
    folder = None
    owned = False
    status = {'year': year, 'pollutant': pollutant, 'source_id': source}
    try:
        row = source_row(inventory_path(year, pollutant), source)
        daily, annual, radius, missing = emissions_for(row, year)
        if annual == 0:
            return dict(status, state='skipped_zero_emissions')
        folder = Path(OUTPUT_BASE_DIR_TEMPLATE.format(year=year, pollutant=pollutant, well_id=source, wellpad_id=source, facility_id=source)).resolve()
        status['folder'] = str(folder)
        fingerprint = task_fingerprint(task)
        if folder.exists():
            done = folder / 'AERMOD_DONE.json'
            archived = None if done.is_file() else consolidated_done(folder)
            if done.is_file() or archived is not None:
                old = json.loads(done.read_text()) if done.is_file() else archived
                previous = dict(old['fingerprint'])
                current = dict(fingerprint)
                old_hash = previous.pop('script_sha256', None)
                current_hash = current.pop('script_sha256', None)
                require(previous == current and old_hash == current_hash, 'Completed case settings changed; use a new output base')
                for info in old['output_files']:
                    require(file_record(info['path']) == info, 'Completed output changed')
                for info in list(old['met_context']['files'].values()) + old['receptor_input_files']:
                    require(file_record(info['path']) == info, 'Completed case meteorology/receptors changed')
                csv_records = [info for info in old['output_files'] if Path(info['path']).suffix.lower() == '.csv']
                require(len(csv_records) == 1, 'Expected one completed daily CSV')
                validate_saved_daily_csv(csv_records[0]['path'], year, old['unique_receptors'])
                removed = []
                for info in old['output_files']:
                    raw = Path(info['path'])
                    if raw.suffix.lower() == '.out':
                        require(raw.resolve().parent == folder and (not raw.is_symlink()), 'Unexpected raw output path')
                        raw.unlink()
                        removed.append(info)
                if removed or old['fingerprint'] != fingerprint:
                    old['fingerprint'] = fingerprint
                    old['output_files'] = csv_records
                    old.setdefault('deleted_raw_outputs', []).extend(removed)
                    old['notes'] = 'Raw OUT deleted after saved CSV validation. Model checks retained in this manifest.'
                    dump(done, old)
                    dump(folder / 'AERMOD_STATUS.json', old)
                compact_success_diagnostics(folder)
                return dict(status, state='skipped_verified_complete')
            raise ValueError('Case folder already exists without AERMOD_DONE.json; preserve it and use a new output base. Another job may own it.')
        folder.parent.mkdir(parents=True, exist_ok=True)
        folder.mkdir()
        owned = True
        status.update(state='started', fingerprint=fingerprint)
        dump(folder / 'AERMOD_STATUS.json', status)
        lat, lon, elevation = coordinates_and_elevation(source, row)
        dist = distances(lat, lon)
        gid = str(grids().iloc[int(np.argmin(dist))]['grid_id'])
        met = meteorology(folder, gid, year)
        receptor_file, points, receptor_inputs = receptors_for(folder, lat, lon, radius)
        stack = stack_for(source)
        stack.update(dict(zip(['easting', 'northing'], utm.from_latlon(lat, lon, force_zone_number=18)[:2])))
        hourly = write_hourly(folder, source, year, daily, stack)
        if WORKFLOW in {'activity', 'production'}:
            inp = generate_aermod_inp(source, year, hourly, str(folder), lat, lon, elevation, pollutant, receptor_file)
        else:
            inp = generate_aermod_inp(source, year, stack, hourly, str(folder), str(folder / 'output_stage2.SFC'), str(folder / 'output_stage2.PFL'), receptor_file, elevation, pollutant)
        require(inp, 'Input writer failed')
        output, model = model_run(inp, folder, year)
        result = daily_csv(output, year, points)
        validate_saved_daily_csv(result, year, len(points))
        raw_record = file_record(output)
        status.update(state='completed', model=model, annual_emissions_tonnes=annual, missing_emission_entries_treated_as_zero=missing, stack_parameters=stack, met_context=met, receptor_input_files=receptor_inputs, unique_receptors=len(points), output_files=[file_record(result)], deleted_raw_outputs=[raw_record], notes='Raw OUT deleted after saved CSV validation. Missing meteorological hours affect daily means. No missing receptor values replaced with zero.')
        require(output.resolve().parent == folder and (not output.is_symlink()), 'Unexpected raw output path')
        output.unlink()
        dump(folder / 'AERMOD_DONE.json', status)
        dump(folder / 'AERMOD_STATUS.json', status)
        compact_success_diagnostics(folder)
        return {**{k: status[k] for k in ['year', 'pollutant', 'source_id', 'folder']}, 'state': 'completed', 'missing_met_hours': model['missing_met_hours_reported']}
    except Exception as exc:
        status.update(state='failed', error=str(exc))
        if folder is not None and owned:
            dump(folder / 'AERMOD_STATUS.json', status)
            (folder / 'failure.txt').write_text(traceback.format_exc())
        return {k: v for k, v in status.items() if k in ['year', 'pollutant', 'source_id', 'folder', 'state', 'error']}

def tasks(log):
    for year in YEARS:
        for pollutant in POLLUTANTS:
            path = Path(inventory_path(year, pollutant))
            if not path.is_file():
                log.write(json.dumps({'year': year, 'pollutant': pollutant, 'state': 'inventory_absent', 'path': str(path)}) + '\n')
                log.flush()
                continue
            df = table(str(path))
            key = id_column()
            require(key in df.columns, f'{path}: missing {key}')
            ids = df[key].map(identifier)
            require(not ids.duplicated().any(), f'Duplicate source IDs in {path}')
            for source in ids:
                if SOURCE_IDS is None or source in {identifier(v) for v in SOURCE_IDS}:
                    yield (int(year), pollutant, source)

def main():
    require(Path(AERMOD_EXECUTABLE).is_file(), 'AERMOD executable not found')
    require(MAX_WORKERS <= int(os.environ.get('SLURM_CPUS_PER_TASK', str(MAX_WORKERS))), 'Workers exceed job CPU allocation')
    require(MAX_WORKERS >= 1 and YEARS and (len(set(YEARS)) == len(YEARS)), 'Check years/workers')
    require(len(set(POLLUTANTS)) == len(POLLUTANTS), 'Duplicate pollutants')
    logdir = Path(LOG_DIRECTORY)
    logdir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '_' + str(os.getpid())
    path = logdir / f'{WORKFLOW}_{run_id}.jsonl'
    failed = 0
    submitted = 0
    with path.open('x') as log, ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        iterator = iter(tasks(log))
        pending = set()
        owners = {}
        finished = False
        while pending or not finished:
            while len(pending) < MAX_WORKERS * 2 and (not finished):
                try:
                    task = next(iterator)
                except StopIteration:
                    finished = True
                    break
                future = pool.submit(run_case, task)
                pending.add(future)
                owners[future] = task
                submitted += 1
            if not pending:
                break
            ready, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in ready:
                task = owners.pop(future)
                try:
                    result = future.result()
                except Exception as exc:
                    result = {'year': task[0], 'pollutant': task[1], 'source_id': task[2], 'state': 'failed', 'error': str(exc)}
                if result.get('state') in {'completed', 'skipped_verified_complete'} and result.get('folder'):
                    try:
                        consolidate_case_audit(Path(result['folder']))
                    except Exception as exc:
                        result['audit_consolidation_error'] = str(exc)
                        print('Audit consolidation incomplete; remaining loose files retained: ' + str(exc), flush=True)
                failed += result['state'] == 'failed'
                log.write(json.dumps(result) + '\n')
                log.flush()
                print(json.dumps(result), flush=True)
    print(f'Submitted {submitted}; failed {failed}; summary: {path}', flush=True)
    require(submitted > 0, 'No cases submitted; inspect inventory paths and selections')
    if failed:
        raise RuntimeError(f'{failed} cases failed; see {path}')
import source_parameters as _sources
from source_parameters import run_case, write_input as generate_aermod_inp, meteorology, receptors as receptors_for, coordinates as coordinates_and_elevation, stack_values as stack_for, fixed_stack as get_stack_params, fingerprint as task_fingerprint, write_hourly, inventory as inventory_path, id_column, emissions_for
_sources.configure(globals())
if __name__ == '__main__':
    main()
