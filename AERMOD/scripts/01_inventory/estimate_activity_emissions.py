"""Calculate emissions from activity and emission factors."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
import argparse
import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path
import sys

def number(value, label):
    try:
        result = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError(f'{label}: invalid number {value!r}') from exc
    if not result.is_finite() or result < 0:
        raise ValueError(f'{label}: must be finite and nonnegative')
    return result

def read_table(filename, fields):
    with Path(filename).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream)
        missing = set(fields) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'{filename}: missing columns {sorted(missing)}')
        rows = list(reader)
    for n, row in enumerate(rows, 2):
        if any((not (row.get(key) or '').strip() for key in fields)):
            raise ValueError(f'{filename}: blank required value on row {n}')
    return rows

def calculate(activity, factors):
    """One row per activity record and factor; unmatched activity is an error."""
    index = {}
    for row in factors:
        key = tuple((row[x] for x in ('component', 'technology', 'activity_unit')))
        bucket = index.setdefault(key, {})
        if row['pollutant'] in bucket:
            raise ValueError(f"Duplicate factor for {key}, {row['pollutant']}")
        number(row['factor_g_per_unit'], 'factor')
        bucket[row['pollutant']] = row
    output = []
    seen = set()
    for row in activity:
        if row['record_id'] in seen:
            raise ValueError(f"Duplicate activity record_id {row['record_id']}")
        seen.add(row['record_id'])
        key = tuple((row[x] for x in ('component', 'technology', 'activity_unit')))
        if key not in index:
            raise ValueError(f"No unit-matched factor for activity {row['record_id']}: {key}")
        amount = number(row['activity_value'], 'activity')
        for pollutant, factor in index[key].items():
            grams = amount * number(factor['factor_g_per_unit'], 'factor')
            output.append({**row, 'pollutant': pollutant, 'emission_g': str(grams), 'emission_kg': str(grams / 1000), 'factor_reference': factor['reference']})
    return output

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--activity', required=True)
    parser.add_argument('--factors', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    activity = read_table(args.activity, ['record_id', 'source_id', 'year', 'component', 'technology', 'activity_unit', 'activity_value'])
    factors = read_table(args.factors, ['component', 'technology', 'activity_unit', 'pollutant', 'factor_g_per_unit', 'reference'])
    output = calculate(activity, factors)
    if not output:
        raise ValueError('No emissions produced; activity or factors are empty')
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
if __name__ == '__main__':
    main()
