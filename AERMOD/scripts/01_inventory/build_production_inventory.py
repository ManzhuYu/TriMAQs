"""Combine leakage and heater emission inventories by well pad and month."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
from pathlib import Path
from collections import defaultdict
from decimal import Decimal
import csv, calendar, hashlib, json, re, zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
BASE = Path(_configured_path('build_production_inventory.base_1'))
LEAK_FOLDER = BASE / 'Production_emission_wellpad2'
HEATER_FILE = BASE / 'emission_result_Heater_2.xlsx'
MAPPING_FILE = BASE / 'production_well_to_pad_mapping_reconciled.csv'
OUTPUT_FOLDER = BASE / 'Production_emission_wellpad3'
YEARS = _setting('build_production_inventory.years')
DROP_EXACT_DUPLICATE_HEATER_ROWS = True
ALLOW_SMALL_UNMATCHED_EXCLUSION = True
MAX_UNMATCHED_FRACTION = 0.0001
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

def require(ok, message):
    if not ok:
        raise ValueError(message)

def key(value):
    s = str(value).strip()
    if not s or s.lower() == 'nan':
        return ''
    try:
        n = Decimal(s)
        if n.is_finite() and n == n.to_integral_value():
            return str(int(n))
    except Exception:
        pass
    return re.sub('[^0-9]', '', s)

def mapping():
    d = pd.read_csv(MAPPING_FILE, dtype=str)
    required = {'well_key', 'pad_id', 'latitude', 'longitude'}
    require(required <= set(d.columns), 'Mapping CSV columns are incorrect')
    lookup = defaultdict(set)
    coords = defaultdict(set)
    for r in d.to_dict('records'):
        w = key(r['well_key'])
        p = key(r['pad_id'])
        require(w and p, 'Blank mapping identifier')
        lookup[w].add(p)
        coords[p].add((float(r['latitude']), float(r['longitude'])))
    return (lookup, {p: tuple(np.mean(list(v), axis=0)) for p, v in coords.items()})

def build():
    require(not OUTPUT_FOLDER.exists(), 'Output already exists; choose a new OUTPUT_FOLDER, never overwrite')
    require(HEATER_FILE.is_file() and MAPPING_FILE.is_file(), 'Copy heater workbook and mapping CSV to configured paths')
    lookup, coords = mapping()
    OUTPUT_FOLDER.mkdir(parents=True)
    totals = []
    unmatched = defaultdict(lambda: {'rows': 0, 'tonnes': 0.0})
    heat = {}
    header_info = {}
    with zipfile.ZipFile(HEATER_FILE) as z:
        strings = [''.join(n.itertext()) for n in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        wb = ET.fromstring(z.read('xl/workbook.xml'))
        rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheets = {s.get('name'): rels[s.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')] for s in wb.find(NS + 'sheets')}
        for pollutant in ['NOx', 'PM2.5', 'CH4', 'VOC']:
            require(pollutant in sheets, f'Missing heater sheet: {pollutant}')
            target = sheets[pollutant]
            target = target.lstrip('/') if target.startswith('/') else 'xl/' + target
            selected = {}
            headers = {}
            seen = set()
            duplicates = 0
            rows = 0
            raw = defaultdict(float)
            dedup = defaultdict(float)
            mapped = defaultdict(float)
            with z.open(target) as stream:
                for event, row in ET.iterparse(stream, events=['end']):
                    if row.tag != NS + 'row':
                        continue
                    values = {}
                    for c in row:
                        v = c.findtext(NS + 'v')
                        col = re.sub('\\d', '', c.get('r'))
                        values[col] = strings[int(v)] if c.get('t') == 's' and v is not None else v
                    if not headers:
                        headers = values
                        permit = next((c for c, h in headers.items() if h == 'PERMIT_NUM'), None)
                        require(permit is not None, 'Heater sheet has no PERMIT_NUM')
                        for c, h in headers.items():
                            m = re.fullmatch('(\\d{2})_\\d{2}_(20\\d{2})_Gas', h or '')
                            if m and int(m[2]) in YEARS:
                                selected[c] = (int(m[2]), int(m[1]))
                        require(selected, 'No requested years in heater sheet')
                        header_info[pollutant] = list(selected.values())
                        row.clear()
                        continue
                    rows += 1
                    vals = {}
                    for c, (year, month) in selected.items():
                        v = values.get(c)
                        n = 0.0 if v is None or v == '' else float(v)
                        require(np.isfinite(n) and n >= 0, 'Invalid heater value')
                        vals[year, month] = n
                        raw[year] += n
                    signature = hashlib.sha256(repr(tuple(sorted(values.items()))).encode()).digest()
                    duplicate = signature in seen
                    seen.add(signature)
                    if duplicate:
                        duplicates += 1
                    if duplicate and DROP_EXACT_DUPLICATE_HEATER_ROWS:
                        row.clear()
                        continue
                    for (year, month), v in vals.items():
                        dedup[year] += v
                    w = key(values.get(permit, ''))
                    pads = lookup.get(w, set())
                    if len(pads) != 1:
                        amount = sum(vals.values())
                        if amount > 0:
                            u = unmatched[pollutant, w, 'unmatched' if not pads else 'ambiguous']
                            u['rows'] += 1
                            u['tonnes'] += amount
                        row.clear()
                        continue
                    pad = next(iter(pads))
                    for (year, month), v in vals.items():
                        if v:
                            bucket = heat.setdefault((pollutant, year, pad), np.zeros(12))
                            bucket[month - 1] += v
                            mapped[year] += v
                    row.clear()
            for year in YEARS:
                totals.append({'pollutant': pollutant, 'year': year, 'raw_workbook_tonnes': raw[year], 'after_duplicate_policy_tonnes': dedup[year], 'mapped_heater_tonnes': mapped[year]})
            print(f'{pollutant}: rows={rows}; exact duplicates={duplicates}; indexed', flush=True)
    pd.DataFrame(totals).to_csv(OUTPUT_FOLDER / 'heater_mass_balance.csv', index=False)
    missing = [dict(pollutant=k[0], well_key=k[1], reason=k[2], **v) for k, v in unmatched.items()]
    pd.DataFrame(missing, columns=['pollutant', 'well_key', 'reason', 'rows', 'tonnes']).to_csv(OUTPUT_FOLDER / 'unmatched_heater_ids.csv', index=False)
    excluded = {p: sum((v['tonnes'] for k, v in unmatched.items() if k[0] == p)) for p in ['NOx', 'PM2.5', 'CH4', 'VOC']}
    denominators = {p: sum((r['after_duplicate_policy_tonnes'] for r in totals if r['pollutant'] == p)) for p in excluded}
    fractions = {p: excluded[p] / denominators[p] if denominators[p] else 0 if not excluded[p] else 1 for p in excluded}
    exclusion_ok = ALLOW_SMALL_UNMATCHED_EXCLUSION and all((v <= MAX_UNMATCHED_FRACTION for v in fractions.values()))
    audit = {'heater_file': str(HEATER_FILE), 'mapping_file': str(MAPPING_FILE), 'units_assumption': 'metric tonnes/month, as supplied; original calculation not independently validated', 'duplicate_policy': 'drop exact rows' if DROP_EXACT_DUPLICATE_HEATER_ROWS else 'retain all rows', 'unmatched_positive_id_groups': len(missing), 'SO2': 'not estimated; no file or synthetic zero inventory generated', 'excluded_heater_tonnes_all_years': excluded, 'excluded_fractions': fractions, 'exclusion_policy': 'Retain unmatched records in audit and exclude only within explicit aggregate mass limit; not physically zero.', 'state': 'blocked_unmatched_ids' if missing and (not exclusion_ok) else 'building'}
    (OUTPUT_FOLDER / 'BUILD_SUMMARY.json').write_text(json.dumps(audit, indent=2))
    require(not missing or exclusion_ok, 'Unmatched emissions exceed permitted exclusion policy. Share unmatched_heater_ids.csv; combined files were NOT created.')
    balances = []
    for year in YEARS:
        leak = {}
        for pol in ['CH4', 'VOC']:
            path = LEAK_FOLDER / f'Emission_{pol}_{year}_WellPad.csv'
            require(path.is_file(), f'Missing confirmed leakage inventory: {path}')
            d = pd.read_csv(path, dtype={'WellPadID': str})
            d['WellPadID'] = d['WellPadID'].map(key)
            require(not d['WellPadID'].duplicated().any(), 'Duplicate pad IDs in leakage inventory')
            records = {}
            for r in d.to_dict('records'):
                vals = np.array([float(r.get(f'{calendar.month_abbr[m]}_{str(year)[-2:]}_Production', 0) or 0) for m in range(1, 13)])
                require(np.isfinite(vals).all() and (vals >= 0).all(), 'Invalid leakage value')
                pad = r['WellPadID']
                records[pad] = vals
                coords[pad] = (float(r['LATITUDE_DECIMAL']), float(r['LONGITUDE_DECIMAL']))
            leak[pol] = records
        for pol in ['NOx', 'PM2.5', 'CH4', 'VOC']:
            heating = {pad: v for (p, y, pad), v in heat.items() if p == pol and y == year}
            leakage = leak.get(pol, {})
            ids = sorted(set(heating) | set(leakage))
            records = []
            for pad in ids:
                v = heating.get(pad, np.zeros(12)) + leakage.get(pad, np.zeros(12))
                lat, lon = coords[pad]
                r = {'WellPadID': pad, 'LATITUDE_DECIMAL': lat, 'LONGITUDE_DECIMAL': lon}
                r.update({f'{calendar.month_abbr[m]}_{str(year)[-2:]}_Production': v[m - 1] for m in range(1, 13)})
                records.append(r)
            out = pd.DataFrame(records)
            columns = ['WellPadID'] + [f'{calendar.month_abbr[m]}_{str(year)[-2:]}_Production' for m in range(1, 13)] + ['LATITUDE_DECIMAL', 'LONGITUDE_DECIMAL']
            out = out.reindex(columns=columns)
            out.to_csv(OUTPUT_FOLDER / f'Emission_{pol}_{year}_WellPad.csv', index=False)
            h = sum((float(v.sum()) for v in heating.values()))
            l = sum((float(v.sum()) for v in leakage.values()))
            total = float(out[columns[1:13]].sum().sum())
            require(abs(total - h - l) <= 1e-08 * max(1, total), 'Combined emission mass did not reconcile')
            balances.append({'year': year, 'pollutant': pol, 'pads': len(ids), 'heater_tonnes': h, 'leak_tonnes': l, 'combined_tonnes': total})
    pd.DataFrame(balances).to_csv(OUTPUT_FOLDER / 'combined_mass_balance.csv', index=False)
    audit.update(state='completed_with_documented_exclusions' if missing else 'completed', output_files=len(YEARS) * 4, years=YEARS)
    (OUTPUT_FOLDER / 'BUILD_SUMMARY.json').write_text(json.dumps(audit, indent=2))
    print('Share BUILD_SUMMARY.json and combined_mass_balance.csv before the AERMOD run.', flush=True)
if __name__ == '__main__':
    build()
