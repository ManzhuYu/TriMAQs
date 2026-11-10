"""Estimate compressor emissions from production and emission factors."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
import os
import re
import pandas as pd
target_years = _setting('estimate_compressor_emissions.years')
pollutants = _setting('estimate_compressor_emissions.pollutants')
compressor_prod_folder = _configured_path('estimate_compressor_emissions.production_folder')
emission_factor_file = _configured_path('estimate_compressor_emissions.emission_factors')
output_base_folder = _configured_path('estimate_compressor_emissions.output_folder')
facility_file = _configured_path('estimate_compressor_emissions.facilities')
for pollutant in pollutants:
    os.makedirs(os.path.join(output_base_folder, pollutant), exist_ok=True)
conv_factor = 1
expected_months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

def normalize_columns(df):
    rename_dict = {}
    for col in df.columns:
        match = re.match('([A-Za-z]+)_\\d+_Production$', col, re.IGNORECASE)
        if match:
            rename_dict[col] = match.group(1).capitalize()
        match2 = re.match('([A-Za-z]+)_\\d+_GAS_QUANTITY$', col, re.IGNORECASE)
        if match2:
            rename_dict[col] = match2.group(1).capitalize()
    return df.rename(columns=rename_dict)
ef_df = pd.read_csv(emission_factor_file)
facility_df = pd.read_csv(facility_file)
facility_df.columns = [col.strip().lower() for col in facility_df.columns]
facility_df['facility_id'] = facility_df['facility_id'].astype(str).str.strip()
for pollutant in pollutants:
    pollutant_folder = os.path.join(output_base_folder, pollutant)
    os.makedirs(pollutant_folder, exist_ok=True)
    for year in target_years:
        print(f'\nProcessing {pollutant} for year {year}...')
        prod_file = os.path.join(compressor_prod_folder, f'compressor_wise_production_{year}.csv')
        if not os.path.exists(prod_file):
            print(f'  Skipping {year} – Production file missing.')
            continue
        prod_df = pd.read_csv(prod_file)
        if 'facility_id' not in prod_df.columns and 'facility_id_compressors' in prod_df.columns:
            prod_df.rename(columns={'facility_id_compressors': 'facility_id'}, inplace=True)
        if 'facility_id' not in prod_df.columns:
            print(f'  Skipping {year} – Missing facility identifier.')
            continue
        prod_df = normalize_columns(prod_df)
        for month in expected_months:
            if month not in prod_df.columns:
                prod_df[month] = 0.0
            else:
                prod_df[month] = pd.to_numeric(prod_df[month], errors='coerce').fillna(0)
        for month in expected_months:
            prod_df[month] = prod_df[month] / 1000000.0
        ef_row = ef_df[ef_df['Year'] == year]
        if ef_row.empty or pollutant not in ef_row.columns:
            print(f'  Emission factor not available for {pollutant} in {year}, skipping.')
            continue
        factor = float(ef_row[pollutant].iloc[0])
        for month in expected_months:
            prod_df[month] = prod_df[month] * factor * conv_factor
        prod_df['facility_id'] = prod_df['facility_id'].astype(str).str.strip()
        merged_df = pd.merge(prod_df, facility_df[['facility_id', 'latitude', 'longitude']], on='facility_id', how='left')
        final_cols = ['facility_id', 'latitude', 'longitude'] + expected_months
        merged_df = merged_df[final_cols]
        output_file = os.path.join(pollutant_folder, f'ONG_emission_{year}.csv')
        merged_df.to_csv(output_file, index=False)
        print(f'  Saved emissions for {pollutant} in {year} to: {output_file}')
