"""Interpolate facility emission inventories between available reporting years."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
import os
import re
import pandas as pd
latest_year = 2023
pollutants = _setting('interpolate_facility_emissions.pollutants')
target_years = _setting('interpolate_facility_emissions.years')
raw_emission_folder = _configured_path('interpolate_facility_emissions.emission_folder')
facility_prod_folder = _configured_path('interpolate_facility_emissions.production_folder')
output_base_folder = _configured_path('interpolate_facility_emissions.output_folder')
for pollutant in pollutants:
    os.makedirs(os.path.join(output_base_folder, pollutant), exist_ok=True)
month_days = {'Jan': 31, 'Feb': 28, 'Mar': 31, 'Apr': 30, 'May': 31, 'Jun': 30, 'Jul': 31, 'Aug': 31, 'Sep': 30, 'Oct': 31, 'Nov': 30, 'Dec': 31}

def get_days_in_month(month, year):
    return 29 if month == 'Feb' and year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else month_days[month]
conv_factor = 0.90718474
expected_months = list(month_days.keys())
emission_rate_dict = {}
for pollutant in pollutants:
    print(f'\nCalculating emission rates for {pollutant}...')
    raw_file = os.path.join(raw_emission_folder, f'{pollutant}_raw_emission.csv')
    if not os.path.exists(raw_file):
        continue
    raw_df = pd.read_csv(raw_file)
    raw_df['facility_id'] = raw_df['facility_id'].astype(str).str.strip()
    raw_em_cols = [col for col in raw_df.columns if re.match('\\d{4}_emission', col)]
    available_years = sorted([int(col.split('_')[0]) for col in raw_em_cols])
    rate_rows = []
    for _, row in raw_df.iterrows():
        facility = row['facility_id']
        year_em_list = [(y, row.get(f'{y}_emission', 0)) for y in available_years if row.get(f'{y}_emission', 0) > 0]
        year_em_list = sorted(year_em_list, key=lambda x: x[0], reverse=True)[:3]
        raw_sum = sum((em for _, em in year_em_list))
        prod_sum = 0.0
        for y, _ in year_em_list:
            prod_file = os.path.join(facility_prod_folder, f'facility_production_{pollutant}_{y}.csv')
            if not os.path.exists(prod_file):
                continue
            prod_df = pd.read_csv(prod_file)
            prod_df['facility_id'] = prod_df['facility_id'].astype(str).str.strip()
            fac_prod = prod_df[prod_df['facility_id'] == facility]
            if fac_prod.empty:
                continue

            def normalize_columns(df):
                rename_dict = {}
                for col in df.columns:
                    match = re.match('([A-Za-z]+)_?\\d{2}_?(production|GAS_QUANTITY)', col, re.IGNORECASE)
                    if match:
                        rename_dict[col] = match.group(1).capitalize()
                return df.rename(columns=rename_dict)
            fac_prod = normalize_columns(fac_prod)
            prod_total = fac_prod[expected_months].sum(axis=1).iloc[0] if all((m in fac_prod.columns for m in expected_months)) else 0.0
            prod_sum += prod_total
        rate = raw_sum / prod_sum if prod_sum > 0 else 0.0
        rate_rows.append({'facility_id': facility, 'emission_rate': rate})
    emission_rate_dict[pollutant] = pd.DataFrame(rate_rows)
common_facilities_dict = {}
for pollutant in pollutants:
    common_ids = set()
    for y in range(2012, 2024):
        prod_file = os.path.join(facility_prod_folder, f'facility_production_{pollutant}_{y}.csv')
        if os.path.exists(prod_file):
            df = pd.read_csv(prod_file)
            df['facility_id'] = df['facility_id'].astype(str).str.strip()
            common_ids.update(df['facility_id'].unique())
    common_facilities_dict[pollutant] = common_ids
for pollutant in pollutants:
    if pollutant not in emission_rate_dict:
        continue
    rate_df = emission_rate_dict[pollutant]
    common_facilities = common_facilities_dict.get(pollutant, set())
    for year in target_years:
        print(f'Interpolating {pollutant} for year {year}...')
        prod_file = os.path.join(facility_prod_folder, f'facility_production_{pollutant}_{year}.csv')
        if not os.path.exists(prod_file):
            continue
        prod_df = pd.read_csv(prod_file)
        prod_df['facility_id'] = prod_df['facility_id'].astype(str).str.strip()

        def normalize_columns(df):
            rename_dict = {}
            for col in df.columns:
                match = re.match('([A-Za-z]+)_?\\d{2}_?(production|GAS_QUANTITY)', col, re.IGNORECASE)
                if match:
                    rename_dict[col] = match.group(1).capitalize()
            return df.rename(columns=rename_dict)
        prod_df = normalize_columns(prod_df)
        missing_months = [m for m in expected_months if m not in prod_df.columns]
        if missing_months:
            continue
        prod_df = prod_df.melt(id_vars=['facility_id'], value_vars=expected_months, var_name='month', value_name='production')
        prod_df['year'] = year
        if year < 2012:
            prod_df = prod_df[prod_df['facility_id'].isin(common_facilities)]
        merged = pd.merge(prod_df, rate_df, on='facility_id', how='inner')
        merged['monthly_emission'] = merged['production'] * merged['emission_rate']
        merged['metric_ton_per_day'] = merged.apply(lambda row: row['monthly_emission'] * conv_factor / get_days_in_month(row['month'], row['year']), axis=1)
        final = merged.groupby(['facility_id', 'month'], as_index=False).agg(metric_ton_per_day=('metric_ton_per_day', 'sum'))
        final = final.pivot(index='facility_id', columns='month', values='metric_ton_per_day').reset_index()
        final.fillna(0, inplace=True)
        if {'latitude', 'longitude'}.issubset(prod_df.columns):
            loc_df = prod_df.drop_duplicates('facility_id')[['facility_id', 'latitude', 'longitude']]
            final = pd.merge(final, loc_df, on='facility_id', how='left')
            final = final[['facility_id', 'latitude', 'longitude'] + expected_months]
        else:
            final = final[['facility_id'] + expected_months]
        output_file = os.path.join(output_base_folder, pollutant, f'ONG_emission_{year}.csv')
        final.to_csv(output_file, index=False)
        print(f'Saved: {output_file}')
