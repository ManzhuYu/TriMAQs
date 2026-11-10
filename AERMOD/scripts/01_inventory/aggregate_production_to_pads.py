"""Aggregate well production-related emissions to well pads."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
import os
import pandas as pd
BASE_IN = _configured_path('aggregate_production_to_pads.well_emission_folder')
BASE_OUT = _configured_path('aggregate_production_to_pads.output_folder')
os.makedirs(BASE_OUT, exist_ok=True)

def convert_production_to_emissions(df, year):
    leak_fraction = 0.003
    kg_per_mcf = 19.3
    voc_to_ch4_ratio = 0.278
    yy = str(year)[-2:]
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    if 'Latitude' in df.columns:
        df['LATITUDE_DECIMAL'] = df['Latitude']
    if 'Longitude' in df.columns:
        df['LONGITUDE_DECIMAL'] = df['Longitude']
    for m in months:
        prod_col = f'{m}_{year}'
        out_col = f'{m}_{yy}_Production'
        if prod_col not in df.columns:
            print(f'Column missing: {prod_col}')
            df[prod_col] = 0
        df[f'{m}_{yy}_CH4'] = df[prod_col] * leak_fraction * kg_per_mcf / 1000
        df[f'{m}_{yy}_VOC'] = df[f'{m}_{yy}_CH4'] * voc_to_ch4_ratio
        df[out_col] = df[f'{m}_{yy}_CH4']
    return df

def create_year_files(year):
    in_file = os.path.join(BASE_IN, f'Wellpad_production_{year}.csv')
    if not os.path.exists(in_file):
        print(f'Missing: {in_file}')
        return
    df = pd.read_csv(in_file)
    df = convert_production_to_emissions(df, year)
    yy = str(year)[-2:]
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    ch4_cols = ['WellPadID'] + [f'{m}_{yy}_Production' for m in months]
    if 'LATITUDE_DECIMAL' in df.columns:
        ch4_cols += ['LATITUDE_DECIMAL', 'LONGITUDE_DECIMAL']
    df[ch4_cols].to_csv(os.path.join(BASE_OUT, f'Emission_CH4_{year}_WellPad.csv'), index=False)
    df_voc = df.copy()
    for m in months:
        out_col = f'{m}_{yy}_Production'
        df_voc[out_col] = df_voc[f'{m}_{yy}_VOC']
    df_voc[ch4_cols].to_csv(os.path.join(BASE_OUT, f'Emission_VOC_{year}_WellPad.csv'), index=False)
    print(f'Saved year {year}: CH4 & VOC files with OLD column format')
for year in range(2004, 2024):
    create_year_files(year)
