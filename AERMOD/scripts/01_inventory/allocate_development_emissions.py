"""Allocate drilling and completion emissions to daily activity intervals."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
import os
import pandas as pd
from datetime import datetime, timedelta
import re
drill_file = _configured_path('allocate_development_emissions.drilling_workbook')
frac_file = _configured_path('allocate_development_emissions.completion_workbook')
WELL_OUTPUT = _configured_path('allocate_development_emissions.well_output')
PAD_OUTPUT = _configured_path('allocate_development_emissions.pad_output')
truck_file = _configured_path('allocate_development_emissions.transport_workbook')
os.makedirs(WELL_OUTPUT, exist_ok=True)
os.makedirs(PAD_OUTPUT, exist_ok=True)
POLLUTANTS = _setting('allocate_development_emissions.pollutants')
YEARS = _setting('allocate_development_emissions.years')
dr = pd.read_excel(drill_file, dtype=str)
fr = pd.read_excel(frac_file, dtype=str)
truck_df = pd.read_excel(truck_file, sheet_name='Results')

def normalize_api14(x):
    s = re.sub('\\D', '', str(x))
    return s[:14].zfill(14) if s else None
dr['Well_ID'] = dr.get('Unformatted_API_UWI_14', dr.get('Unformatted_API_UWI')).apply(normalize_api14)
fr['Well_ID'] = fr.get('Unformatted_API_UWI_14', fr.get('Unformatted_API_UWI')).apply(normalize_api14)
truck_df['Well_ID'] = truck_df['Unformatted_API_UWI_14'].astype(str).str.zfill(14)
dr = dr.merge(truck_df[['Well_ID', 'Construction_NOx', 'Drilling Rig_NOx', 'Construction_PM2.5', 'Drilling Rig_PM2.5', 'Construction_VOC', 'Drilling Rig_VOC']], on='Well_ID', how='left')
fr = fr.merge(truck_df[['Well_ID', 'Completion Rig_NOx', 'Completion Material_NOx', 'Completion Rig_PM2.5', 'Completion Material_PM2.5', 'Completion Rig_VOC', 'Completion Material_VOC']], on='Well_ID', how='left')
num_cols_dr = [c for c in dr.columns if any((x in c for x in ['NOx', 'PM2.5', 'VOC', 'Methane', 'PM2_5']))]
num_cols_fr = [c for c in fr.columns if any((x in c for x in ['NOx', 'PM2.5', 'VOC', 'CH4', 'PM2_5']))]
dr[num_cols_dr] = dr[num_cols_dr].apply(pd.to_numeric, errors='coerce')
fr[num_cols_fr] = fr[num_cols_fr].apply(pd.to_numeric, errors='coerce')
for df, datecol, label in [(dr, 'DrillingDate', 'Drilling'), (fr, 'Job_Date', 'Fracking')]:
    if datecol in df.columns:
        df[datecol] = pd.to_datetime(df[datecol], errors='coerce')
        yr = df[datecol].dt.year
        before = len(df)
        df = df[yr.between(2004, 2023, inclusive='both')].copy()
        after = len(df)
        print(f'{label}: filtered {before - after} rows outside 2004–2023')
        if label == 'Drilling':
            dr = df
        else:
            fr = df
print('=== Drilling Columns ===')
print(dr.columns.tolist()[:15], '…', len(dr.columns))
print('\n=== Completion Columns ===')
print(fr.columns.tolist()[:15], '…', len(fr.columns))
print(f'\nRows (2004–2023): Drilling={len(dr)}, Completion={len(fr)}')
from datetime import datetime, timedelta
import pandas as pd

def spread_total_over_days(total, start, end):
    """
    Spread a total emission mass evenly over each day between start and end (exclusive of end).
    Returns a dictionary: {datetime.date: daily_value}
    """
    if pd.isna(start) or pd.isna(end):
        return {}
    start = pd.to_datetime(start)
    end = pd.to_datetime(end)
    n_days = max((end - start).days, 1)
    if n_days <= 0:
        return {}
    daily_val = float(total) / n_days
    return {start + timedelta(days=i): daily_val for i in range(n_days)}

def clip_to_year(d2v, year):
    """
    Filters a date→value dict to include only keys (dates) within the target year.
    """
    start = datetime(year, 1, 1)
    end = datetime(year, 12, 31)
    return {d: v for d, v in d2v.items() if start.date() <= d.date() <= end.date()}

def build_year_matrix(rows, year):
    """
    Convert list of dicts [{Well_ID, Latitude, Longitude, per_day}] 
    into a DataFrame with 365 (or 366) daily emission columns for the given year.
    """
    days = pd.date_range(f'{year}-01-01', f'{year}-12-31', freq='D')
    out = []
    for r in rows:
        row = {'Well_ID': r['Well_ID'], 'Latitude': r['Latitude'], 'Longitude': r['Longitude']}
        vals = r['per_day']
        for d in days:
            row[str(d.date())] = vals.get(d, 0.0)
        out.append(row)
    return pd.DataFrame(out)

def check_spread_example():
    """Example diagnostic: verify even spread and clipping behavior."""
    test = spread_total_over_days(10, '2020-01-01', '2020-01-06')
    print('Spread example:', test)
    clipped = clip_to_year(test, 2020)
    print('Clipped example:', clipped)
check_spread_example()
from datetime import timedelta
DR_P50 = {'NOx': 'NOx_Emission_DrillingTime_P50', 'SO2': 'SO2_Emission_DrillingTime_P50', 'PM2.5': 'PM2_5_Emission_DrillingTime_P50', 'VOC': 'VOC_Emission_DrillingTime_P50'}
INIT_COL = {'NOx': 'NOx_InitialEmission', 'SO2': 'SO2_InitialEmission', 'PM2.5': 'PM2_5_InitialEmission', 'VOC': 'VOC_InitialEmission'}
FR_HP50 = {'NOx': 'NOx_Emission_HorsepowerOnly_P50', 'SO2': 'SO2_Emission_HorsepowerOnly_P50', 'PM2.5': 'PM2_5_Emission_HorsepowerOnly_P50', 'VOC': 'VOC_Emission_HorsepowerOnly_P50'}
fr_w = fr.set_index('Well_ID', drop=False)
for pollutant in POLLUTANTS:
    well_dir = os.path.join(WELL_OUTPUT, pollutant)
    os.makedirs(well_dir, exist_ok=True)
    perwell_daily = {}
    for _, r in dr.iterrows():
        wid, lat, lon = (r['Well_ID'], r.get('Latitude'), r.get('Longitude'))
        if wid not in perwell_daily:
            perwell_daily[wid] = {'lat': lat, 'lon': lon, 'd2v': {}}
        init_col = INIT_COL.get(pollutant)
        if pd.notna(r.get('DrillingDate')) and pd.notna(r.get('EmissionDate')):
            init_start = pd.to_datetime(r['DrillingDate'])
            init_end = pd.to_datetime(r['EmissionDate'])
            dur = min(max((init_end - init_start).days, 1), 3)
            init_end = init_start + timedelta(days=dur)
            base_init = float(r.get(init_col, 0.0) or 0.0) if init_col else 0.0
            truck_init = 0.0
            if pollutant == 'NOx':
                truck_init = float(r.get('Construction_NOx', 0.0) or 0.0)
            elif pollutant == 'PM2.5':
                truck_init = float(r.get('Construction_PM2.5', 0.0) or 0.0)
            elif pollutant == 'VOC':
                truck_init = float(r.get('Construction_VOC', 0.0) or 0.0)
            total_init = base_init + truck_init
            if total_init > 0:
                init_map = spread_total_over_days(total_init, init_start, init_end)
                for d, v in init_map.items():
                    perwell_daily[wid]['d2v'][d] = perwell_daily[wid]['d2v'].get(d, 0.0) + v
        dr_col = DR_P50.get(pollutant)
        if pd.notna(r.get('EmissionDate')):
            dr_start = pd.to_datetime(r['EmissionDate'])
            dr_days = int(pd.to_numeric(r.get('DrillingTime', 1), errors='coerce') or 1)
            dr_days = max(dr_days, 1)
            dr_end = dr_start + timedelta(days=dr_days)
            base_dr = float(r.get(dr_col, 0.0) or 0.0) if dr_col else 0.0
            truck_dr = 0.0
            if pollutant == 'NOx':
                truck_dr = float(r.get('Drilling Rig_NOx', 0.0) or 0.0)
            elif pollutant == 'PM2.5':
                truck_dr = float(r.get('Drilling Rig_PM2.5', 0.0) or 0.0)
            elif pollutant == 'VOC':
                truck_dr = float(r.get('Drilling Rig_VOC', 0.0) or 0.0)
            total_dr = base_dr + truck_dr
            if total_dr > 0:
                dr_map = spread_total_over_days(total_dr, dr_start, dr_end)
                for d, v in dr_map.items():
                    perwell_daily[wid]['d2v'][d] = perwell_daily[wid]['d2v'].get(d, 0.0) + v
        if pollutant == 'CH4':
            ch4_dr = float(r.get('Methane_Emission', 0.0) or 0.0)
            if ch4_dr > 0 and pd.notna(r.get('DrillingDate')):
                ch4_start = pd.to_datetime(r['DrillingDate'])
                ch4_days = int(pd.to_numeric(r.get('DrillingTime', 1), errors='coerce') or 1)
                ch4_days = max(ch4_days, 1)
                ch4_end = ch4_start + timedelta(days=ch4_days)
                ch4_map = spread_total_over_days(ch4_dr, ch4_start, ch4_end)
                for d, v in ch4_map.items():
                    perwell_daily[wid]['d2v'][d] = perwell_daily[wid]['d2v'].get(d, 0.0) + v
        if wid in fr_w.index:
            fr_subset = fr_w.loc[wid]
            if isinstance(fr_subset, pd.Series):
                fr_subset = fr_subset.to_frame().T
            for _, fr_row in fr_subset.iterrows():
                if pd.isna(fr_row.get('Job_Date')):
                    continue
                fr_start = pd.to_datetime(fr_row['Job_Date'])
                fr_days = int(pd.to_numeric(fr_row.get('CompletionTime_DAYS', 1), errors='coerce') or 1)
                fr_days = max(fr_days, 1)
                fr_end = fr_start + timedelta(days=fr_days)
                fr_col = FR_HP50.get(pollutant)
                base_fr = float(fr_row.get(fr_col, 0.0) or 0.0) if fr_col else 0.0
                truck_fr = 0.0
                if pollutant == 'NOx':
                    truck_fr = float(fr_row.get('Completion Rig_NOx', 0.0) or 0.0) + float(fr_row.get('Completion Material_NOx', 0.0) or 0.0)
                elif pollutant == 'PM2.5':
                    truck_fr = float(fr_row.get('Completion Rig_PM2.5', 0.0) or 0.0) + float(fr_row.get('Completion Material_PM2.5', 0.0) or 0.0)
                elif pollutant == 'VOC':
                    truck_fr = float(fr_row.get('Completion Rig_VOC', 0.0) or 0.0) + float(fr_row.get('Completion Material_VOC', 0.0) or 0.0)
                total_fr = base_fr + truck_fr
                if total_fr > 0:
                    fr_map = spread_total_over_days(total_fr, fr_start, fr_end)
                    for d, v in fr_map.items():
                        perwell_daily[wid]['d2v'][d] = perwell_daily[wid]['d2v'].get(d, 0.0) + v
                if pollutant in ['CH4', 'NOx', 'VOC']:
                    vent_ch4 = float(fr_row.get('CH4_CompletionVentEmission', 0.0) or 0.0)
                    vent_nox = float(fr_row.get('NOx_CompletionVentEmission', 0.0) or 0.0)
                    vent_voc = vent_ch4 / 25.7 if vent_ch4 > 0 else 0.0
                    add_total = 0.0
                    if pollutant == 'CH4':
                        add_total = vent_ch4
                    elif pollutant == 'NOx':
                        add_total = vent_nox
                    elif pollutant == 'VOC':
                        add_total = vent_voc
                    if add_total > 0:
                        vmap = spread_total_over_days(add_total, fr_start, fr_end)
                        for d, v in vmap.items():
                            perwell_daily[wid]['d2v'][d] = perwell_daily[wid]['d2v'].get(d, 0.0) + v
    for year in YEARS:
        rows = []
        for wid, info in perwell_daily.items():
            d2v_year = clip_to_year(info['d2v'], year)
            if not d2v_year:
                continue
            rows.append({'Well_ID': wid, 'Latitude': info['lat'], 'Longitude': info['lon'], 'per_day': d2v_year})
        df = build_year_matrix(rows, year) if rows else build_year_matrix([], year)
        out_path = os.path.join(well_dir, f'{year}_emission.csv')
        df.to_csv(out_path, index=False)
    print(f'✅ Finished {pollutant}')
import re
import numpy as np
import pandas as pd
import os
mapping_file = _configured_path('allocate_development_emissions.pad_mapping')
try:
    map_df = pd.read_csv(mapping_file, dtype=str)
except Exception as e:
    raise RuntimeError(f'❌ Could not read mapping file: {e}')
map_df.columns = [c.strip().lower() for c in map_df.columns]
map_df.rename(columns={'wellid': 'Well_ID', 'well_id': 'Well_ID', 'wellpadid': 'Pad_ID', 'well_pad_id': 'Pad_ID'}, inplace=True)
if 'Well_ID' not in map_df.columns or 'Pad_ID' not in map_df.columns:
    raise KeyError("Mapping file must have columns 'Well_ID' and 'Pad_ID'")
print(f'✅ Loaded mapping file: {len(map_df)} records')
print(map_df.head(3))
for pol in POLLUTANTS:
    src_base = os.path.join(WELL_OUTPUT, pol)
    out_base = os.path.join(PAD_OUTPUT, pol)
    os.makedirs(out_base, exist_ok=True)
    print(f'\n=== Aggregating {pol} to pad level using mapping ===')
    for year in YEARS:
        f = os.path.join(src_base, f'{year}_emission.csv')
        if not os.path.exists(f):
            print(f'⚠️  Missing file for {pol} {year}, skipping…')
            continue
        try:
            df = pd.read_csv(f, dtype={'Well_ID': str})
        except pd.errors.EmptyDataError:
            print(f'⚠️  Empty file (no columns): {pol} {year}, skipping…')
            continue
        except Exception as e:
            print(f'⚠️  Could not read {pol} {year}: {e}')
            continue
        if df.empty or 'Well_ID' not in df.columns:
            print(f'⚠️  No valid data in {pol} {year}, skipping…')
            continue
        df = df.merge(map_df[['Well_ID', 'Pad_ID']], on='Well_ID', how='left')
        if df['Pad_ID'].isna().all():
            print(f'⚠️  No Pad_ID match for {pol} {year}, skipping…')
            continue
        date_cols = [c for c in df.columns if re.match('\\d{4}-\\d{2}-\\d{2}', c)]
        if not date_cols:
            print(f'⚠️  No date columns for {pol} {year}, skipping…')
            continue
        df['Latitude'] = pd.to_numeric(df['Latitude'], errors='coerce')
        df['Longitude'] = pd.to_numeric(df['Longitude'], errors='coerce')
        for c in date_cols:
            df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0.0)
        agg_spec = {dc: 'sum' for dc in date_cols}
        agg_spec.update({'Latitude': 'mean', 'Longitude': 'mean'})
        df_agg = df.groupby('Pad_ID', as_index=False).agg(agg_spec)
        df_agg.rename(columns={'Pad_ID': 'Well_ID'}, inplace=True)
        out_path = os.path.join(out_base, f'{year}_emission.csv')
        df_agg.to_csv(out_path, index=False)
    print(f'✅ Completed pad-level aggregation for {pol}')
