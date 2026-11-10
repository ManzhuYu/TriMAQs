"""Check AERMET results across all selected years and grids."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
from pathlib import Path
from datetime import datetime, timedelta
from collections import Counter
from decimal import Decimal
import calendar
import csv
import json
import re
output_root = Path(_configured_path('check_aermet.output_root_1'))
grid_csv_file = Path(_configured_path('check_aermet.grid_csv_file_2'))
years = _setting('check_aermet.years')
grid_ids = None
verify_hourly_calendar = False

def require(ok, message):
    if not ok:
        raise ValueError(message)

def load_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def normalize_id(value):
    number = Decimal(str(value))
    require(number.is_finite() and number == number.to_integral_value(), f'Invalid grid ID: {value}')
    return str(int(number))

def write_csv(path, fields, rows):
    with path.open('x', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

def hourly_calendar(path, year, surface):
    expected = {datetime(year, 1, 1) + timedelta(hours=h) for h in range(1, (8784 if calendar.isleap(year) else 8760) + 1)}
    seen = set()
    with path.open(errors='strict') as f:
        if surface:
            header = next(f, '')
            require('VERSION:' in header, 'SFC header missing')
        for line in f:
            p = line.split()
            if not p:
                continue
            require(len(p) >= (25 if surface else 11), 'Truncated met record')
            y = int(p[0])
            y = y + 2000 if y < 100 else y
            hour = int(p[4 if surface else 3])
            require(1 <= hour <= 24, f'Invalid hour {hour}')
            key = datetime(y, int(p[1]), int(p[2])) + timedelta(hours=hour)
            if surface:
                require(key not in seen, f'Duplicate SFC hour {key}')
            seen.add(key)
    require(seen == expected, f'Calendar mismatch: missing={len(expected - seen)}, unexpected={len(seen - expected)}')

def inspect_case(folder, year, gid):
    row = {'year': year, 'grid_id': gid, 'status': 'REVIEW', 'folder': str(folder), 'expected_hours': 8784 if calendar.isleap(year) else 8760, 'audited_hours': '', 'noncalm_invalid_hours': '', 'manifest_state': '', 'original_error': '', 'issues': ''}
    issues = []
    if not folder.is_dir():
        row.update(status='MISSING_FOLDER', issues='Normal grid folder absent')
        return row
    required = ['DONE.json', 'manifest.json', 'output_audit.json', 'output_stage2.SFC', 'output_stage2.PFL', 'AERMET_REPORT.RPT', 'stage2.rpt']
    for name in required:
        p = folder / name
        if not p.is_file() or p.stat().st_size == 0:
            issues.append('missing/empty ' + name)
    try:
        manifest = load_json(folder / 'manifest.json') if (folder / 'manifest.json').is_file() else {}
        row['manifest_state'] = manifest.get('state', 'absent')
        row['original_error'] = manifest.get('error', '')
        if manifest.get('state') != 'completed_pending_scientific_review':
            issues.append('manifest not completed')
        if normalize_id(manifest.get('grid_id', gid)) != gid or int(manifest.get('year', year)) != year:
            issues.append('manifest grid/year mismatch')
        done = load_json(folder / 'DONE.json') if (folder / 'DONE.json').is_file() else {}
        if done:
            if done.get('state') != 'completed_pending_scientific_review':
                issues.append('DONE state not completed')
            if normalize_id(done.get('grid_id', '-1')) != gid or int(done.get('year', -1)) != year:
                issues.append('DONE grid/year mismatch')
        audit = load_json(folder / 'output_audit.json') if (folder / 'output_audit.json').is_file() else {}
        row['audited_hours'] = audit.get('hours', '')
        row['noncalm_invalid_hours'] = audit.get('noncalm_invalid_hours', '')
        if audit.get('hours') != row['expected_hours']:
            issues.append('annual audit hour count incorrect/unavailable')
        if audit.get('noncalm_invalid_hours') != 0:
            issues.append('invalid noncalm hours present/unavailable')
        if done.get('output_audit') != audit:
            issues.append('DONE and current audit differ/unavailable')
        provenance = done.get('output_files', {})
        for name in ('output_stage2.SFC', 'output_stage2.PFL'):
            p = folder / name
            saved = provenance.get(name, {})
            if p.is_file():
                stat = p.stat()
                if saved.get('size') != stat.st_size or saved.get('mtime_ns') != stat.st_mtime_ns:
                    issues.append(name + ' differs from completion record/unrecorded')
                if verify_hourly_calendar:
                    try:
                        hourly_calendar(p, year, name.endswith('SFC'))
                    except Exception as exc:
                        issues.append(name + ': ' + str(exc))
    except Exception as exc:
        issues.append('metadata read error: ' + str(exc))
    row['status'] = 'COMPLETE' if not issues else 'REVIEW'
    row['issues'] = '; '.join(issues)
    return row

def main():
    require(output_root.is_dir(), f'AERMET output root not found: {output_root}')
    require(years and len(years) == len(set(years)), 'Empty/duplicate years')
    with grid_csv_file.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        require('grid_ID' in (reader.fieldnames or []), 'Grid CSV needs grid_ID column')
        ids = [normalize_id(r['grid_ID']) for r in reader]
    require(ids and len(ids) == len(set(ids)), 'Empty/duplicate grid IDs in CSV')
    if grid_ids is not None:
        wanted = {normalize_id(g) for g in grid_ids}
        require(wanted and wanted <= set(ids), 'Selected grid IDs not in CSV')
        ids = [g for g in ids if g in wanted]
    expected_ids = set(ids)
    report = output_root / ('CHECK_COMPLETENESS_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    report.mkdir(exist_ok=False)
    all_rows = []
    annual = []
    extras = []
    archives = 0
    for year in years:
        year_rows = [inspect_case(output_root / str(year) / f'Grid_{gid}', year, gid) for gid in ids]
        all_rows.extend(year_rows)
        counts = Counter((r['status'] for r in year_rows))
        annual.append({'year': year, 'expected': len(ids), 'complete': counts['COMPLETE'], 'review': counts['REVIEW'], 'missing_folders': counts['MISSING_FOLDER']})
        directory = output_root / str(year)
        if directory.is_dir():
            for p in directory.iterdir():
                if not p.is_dir():
                    continue
                match = re.fullmatch('Grid_(\\d+)', p.name)
                if match and match[1] not in expected_ids:
                    extras.append({'year': year, 'folder': str(p), 'reason': 'normal folder outside selected expected IDs'})
                elif p.name.startswith('Grid_') and (not match):
                    archives += 1
        print(f"{year}: {counts['COMPLETE']}/{len(ids)} complete; {counts['REVIEW']} review; {counts['MISSING_FOLDER']} missing", flush=True)
    issues = [r for r in all_rows if r['status'] != 'COMPLETE']
    fields = list(all_rows[0])
    write_csv(report / 'all_cases.csv', fields, all_rows)
    write_csv(report / 'issues_only.csv', fields, issues)
    write_csv(report / 'year_summary.csv', list(annual[0]), annual)
    write_csv(report / 'unexpected_folders.csv', ['year', 'folder', 'reason'], extras)
    summary = {'years': years, 'expected_grids_per_year': len(ids), 'expected_cases': len(all_rows), 'complete_cases': len(all_rows) - len(issues), 'review_cases': sum((r['status'] == 'REVIEW' for r in issues)), 'missing_folders': sum((r['status'] == 'MISSING_FOLDER' for r in issues)), 'all_runs_complete': not issues, 'archives_ignored': archives, 'fresh_calendar_scan': verify_hourly_calendar, 'limits': 'Completion inventory using runner audits/provenance; not independent scientific validation. No missing hours waived.', 'report_folder': str(report)}
    with (report / 'summary.json').open('x', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))
    print('Share summary.json and issues_only.csv. Model results were not modified.')
    return summary
if __name__ == '__main__':
    main()
