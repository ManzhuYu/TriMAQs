"""Calculate daily and annual source contributions relative to ML predictions."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting, relative_metadata
from pathlib import Path
from datetime import date, datetime, timedelta
import calendar, json, gc
import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree
BASE = Path(_configured_path('create_source_contributions.base_1'))
INPUT = BASE / 'outputs/concentrations'
OUTPUT = BASE / 'outputs/source_contributions'
ML = BASE / 'inputs/ml_predictions'
PM_ORIGINAL = BASE / 'inputs/ml_predictions/PM25'
YEARS = _setting('create_source_contributions.years')
POLLUTANTS = _setting('create_source_contributions.pollutants')
SOURCES = ['Tier_1', 'Tier_2', 'UOGD']
UNITS = {'NO2': 'ppb', 'SO2': 'ppb', 'PM2.5': 'ug/m3'}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def read(var):
    return np.ma.filled(var[:], np.nan)

def daily_ml_path(pol, stamp):
    if pol == 'PM2.5':
        return PM_ORIGINAL / f'prediction_{stamp}.nc'
    return ML / pol / '1km' / f'{stamp}.nc'

def load_ml(path, study_points, cache, pol):
    """Map using ALL geographic coordinates, independent of missing predictions."""
    with Dataset(path) as ds:
        names = [n for n in ['Prediction', 'prediction'] if n in ds.variables]
        require(len(names) == 1, 'Expected one ML prediction variable: ' + str(path))
        variable = ds.variables[names[0]]
        prediction = np.squeeze(read(variable))
        lat = np.squeeze(read(ds.variables['lat']))
        lon = np.squeeze(read(ds.variables['lon']))
        if lat.ndim == 1 and lon.ndim == 1:
            lon, lat = np.meshgrid(lon, lat)
        require(prediction.ndim == 2 and prediction.shape == lat.shape == lon.shape, 'ML spatial shape mismatch: ' + str(path))
        recorded = str(getattr(variable, 'units', '')).lower().replace('Âµ', 'u').replace('Î¼', 'u').replace(' ', '')
        if recorded:
            allowed = {'ppb', 'ppbv'} if UNITS[pol] == 'ppb' else {'ug/m3', 'ug/m^3', 'ugm-3', 'ugm^-3', 'ug/m**3'}
            require(recorded in allowed, 'ML units differ from required units: ' + str(path) + ' ' + recorded)
        if cache.get('shape') != lat.shape or not np.array_equal(cache.get('lat'), lat, equal_nan=True) or (not np.array_equal(cache.get('lon'), lon, equal_nan=True)):
            valid = np.isfinite(lat) & np.isfinite(lon)
            positions = np.flatnonzero(valid)
            require(len(positions) > 0, 'No finite ML coordinates')
            distance, index = cKDTree(np.column_stack([lat[valid], lon[valid]])).query(study_points)
            require(distance.max() <= 0.0001 and len(np.unique(index)) == len(study_points), 'ML coordinates not aligned one-to-one: ' + str(path))
            cache.update(shape=lat.shape, lat=lat.copy(), lon=lon.copy(), positions=positions[index], maximum_distance=float(distance.max()))
        return prediction.ravel()[cache['positions']].astype(np.float64)

def check_pm_overlap(year, points):
    """Check first/middle/last available old PM files against complete originals."""
    days = 366 if calendar.isleap(year) else 365
    overlap = []
    for day in range(days):
        stamp = (date(year, 1, 1) + timedelta(days=day)).isoformat()
        old = ML / 'PM2.5/1km' / f'{stamp}.nc'
        new = daily_ml_path('PM2.5', stamp)
        if old.is_file() and new.is_file():
            overlap.append((old, new))
    checks = []
    cache1 = {}
    cache2 = {}
    for i in sorted({0, len(overlap) // 2, len(overlap) - 1}) if overlap else []:
        old, new = overlap[i]
        a = load_ml(old, points, cache1, 'PM2.5')
        b = load_ml(new, points, cache2, 'PM2.5')
        require(np.allclose(a, b, rtol=1e-06, atol=1e-08, equal_nan=True), 'PM original differs from older file: ' + str(old) + '; inspect before selecting a model version')
        finite = np.isfinite(a) & np.isfinite(b)
        checks.append(dict(old=str(old), selected=str(new), maximum_difference=float(np.max(abs(a[finite] - b[finite]))) if finite.any() else None))
    return dict(overlapping_dates=len(overlap), sample_checks=checks, limits='First/middle/last overlapping date per year; not proof that every overlapping file is identical.')

def fractions(concentrations, denominator):
    valid = np.isfinite(denominator) & (denominator > 0) & np.isfinite(concentrations).all(axis=0)
    result = np.full(concentrations.shape, np.nan, dtype=np.float64)
    np.divide(concentrations, denominator[None, :], out=result, where=valid[None, :])
    return (result, valid)

def raster(values, inside, shape):
    array = np.full(values.shape[:-1] + shape, np.nan, dtype=np.float32)
    array[..., inside] = values
    return array

def build_year(year):
    final = OUTPUT / f'uogd_source_apportionment_{year}.nc'
    temporary = OUTPUT / f'uogd_source_apportionment_{year}.partial.nc'
    require(not final.exists() and (not temporary.exists()), 'Existing final/partial output; inspect before retry: ' + str(final))
    days = 366 if calendar.isleap(year) else 365
    with Dataset(INPUT / f'aermod_daily_conc_{year}.nc') as inp:
        require(list(inp.variables['source'][:]) == SOURCES, 'Unexpected AERMOD source labels')
        require(len(inp.dimensions['time']) == days, 'Unexpected AERMOD calendar')
        expected_dates = f'days since {year}-01-01 00:00:00'
        require(inp.variables['time'].units == expected_dates and np.array_equal(read(inp.variables['time']), np.arange(days)), 'Unexpected time coordinates')
        mask = read(inp.variables['study_mask']).astype(bool)
        inside = ~mask
        shape = mask.shape
        lat = read(inp.variables['lat'])
        lon = read(inp.variables['lon'])
        points = np.column_stack([lat[inside], lon[inside]])
        overlap = check_pm_overlap(year, points)
        audit = []
        with Dataset(temporary, 'w', format='NETCDF4') as out:
            for key in ['time', 'source', 'y', 'x']:
                out.createDimension(key, len(inp.dimensions[key]))
            for key in ['time', 'source', 'y', 'x', 'lat', 'lon', 'study_mask', 'grid_ID']:
                v = inp.variables[key]
                dtype = str if key == 'source' else v.dtype
                attrs = {k: v.getncattr(k) for k in v.ncattrs() if k != '_FillValue'}
                kwargs = {'fill_value': v.getncattr('_FillValue')} if '_FillValue' in v.ncattrs() else {}
                new = out.createVariable(key, dtype, v.dimensions, **kwargs)
                new[:] = v[:]
                new.setncatts(attrs)
            out.title = f'UOGD contribution relative to ML-predicted total concentrations ({year})'
            out.input_aermod = relative_metadata(str(INPUT / f'aermod_daily_conc_{year}.nc'))
            out.created = datetime.now().isoformat()
            out.ml_units_basis = 'NO2 and SO2: ppb; PM2.5: ug/m3.'
            out.source_definition = 'Tier_1, Tier_2 and UOGD; UOGD is their sum, not an independent category.'
            out.missing_policy = 'Missing dates/cells, nonfinite or nonpositive ML denominator yield NaN; no filling, clipping or interpolation. Coordinate matching only within 1e-4 degree with one-to-one mapping.'
            out.annual_definition = '*_annual_mean_fraction=mean of valid daily ratios. *_annual_ratio=sum of UOGD or tier daily concentrations / sum of ML on the same valid days. Neither imputes missing days.'
            out.crs = getattr(inp, 'crs', '')
            out.parameter_review = getattr(inp, 'assumptions', '')
            for pol in POLLUTANTS:
                name = pol.replace('.', '')
                source_var = inp.variables[name]
                require(source_var.units == UNITS[pol], 'AERMOD concentration units mismatch')
                options = dict(zlib=True, complevel=4, fill_value=np.nan)
                daily = out.createVariable(name + '_fraction', 'f4', ('time', 'source', 'y', 'x'), chunksizes=(1, 1, *shape), **options)
                daily.units = '1'
                daily.coordinates = 'lat lon'
                daily.long_name = f'Daily {pol} source concentration divided by ML total concentration'
                ml = out.createVariable(name + '_ML', 'f4', ('time', 'y', 'x'), chunksizes=(1, *shape), **options)
                ml.units = UNITS[pol]
                ml.coordinates = 'lat lon'
                mean = out.createVariable(name + '_annual_mean_fraction', 'f4', ('source', 'y', 'x'), **options)
                mean.units = '1'
                mean.coordinates = 'lat lon'
                ratio = out.createVariable(name + '_annual_ratio', 'f4', ('source', 'y', 'x'), **options)
                ratio.units = '1'
                ratio.coordinates = 'lat lon'
                count_var = out.createVariable(name + '_valid_days', 'i2', ('y', 'x'), fill_value=-1, zlib=True, complevel=4)
                count_var.long_name = 'Days with finite source concentration and finite positive ML prediction'
                count = np.zeros(len(points), dtype=np.int32)
                sum_fraction = np.zeros((3, len(points)), dtype=float)
                sum_concentration = np.zeros_like(sum_fraction)
                sum_ml = np.zeros(len(points))
                cache = {}
                missing = []
                nonpositive = 0
                missing_cells = 0
                above_one = 0
                for day in range(days):
                    stamp = (date(year, 1, 1) + timedelta(days=day)).isoformat()
                    path = daily_ml_path(pol, stamp)
                    if not path.is_file():
                        missing.append(stamp)
                        continue
                    denominator = load_ml(path, points, cache, pol)
                    conc = np.ma.filled(source_var[day], np.nan)[:, inside].astype(float)
                    require(np.isfinite(conc).all() and (conc >= 0).all(), 'Invalid AERMOD study values')
                    frac, valid = fractions(conc, denominator)
                    require(np.isfinite(frac[:, valid]).all(), 'Invalid finite-denominator fraction')
                    daily[day] = raster(frac, inside, shape)
                    ml[day] = raster(denominator, inside, shape)
                    count += valid
                    sum_fraction[:, valid] += frac[:, valid]
                    sum_concentration[:, valid] += conc[:, valid]
                    sum_ml[valid] += denominator[valid]
                    nonpositive += int((np.isfinite(denominator) & (denominator <= 0)).sum())
                    missing_cells += int((~np.isfinite(denominator)).sum())
                    above_one += int((frac[2] > 1).sum())
                average = np.full_like(sum_fraction, np.nan)
                annual_ratio = np.full_like(sum_fraction, np.nan)
                np.divide(sum_fraction, count[None, :], out=average, where=count[None, :] > 0)
                np.divide(sum_concentration, sum_ml[None, :], out=annual_ratio, where=sum_ml[None, :] > 0)
                mean[:] = raster(average, inside, shape)
                ratio[:] = raster(annual_ratio, inside, shape)
                counts = np.full(shape, -1, dtype=np.int16)
                counts[inside] = count
                count_var[:] = counts
                audit.append(dict(pollutant=pol, available_days=days - len(missing), missing_dates=missing, study_cells_without_valid_days=int((count == 0).sum()), minimum_valid_days=int(count.min()), maximum_valid_days=int(count.max()), nonpositive_cell_days=nonpositive, missing_cell_days_on_available_dates=missing_cells, uogd_fraction_above_one_cell_days=above_one, maximum_coordinate_distance=cache.get('maximum_distance'), ml_template=str(daily_ml_path(pol, 'YYYY-MM-DD'))))
                print(year, pol, 'available:', days - len(missing), 'of', days, flush=True)
                del sum_fraction, sum_concentration, sum_ml
                gc.collect()
            out.coverage_audit = json.dumps(relative_metadata(audit))
            out.pm_overlap_check = json.dumps(relative_metadata(overlap))
    verify(temporary)
    temporary.rename(final)
    return dict(year=year, state='completed', output=str(final), size=final.stat().st_size, coverage=audit, pm_overlap=overlap)

def verify(path):
    """Recompute annual means/counts from saved daily fraction values."""
    with Dataset(path) as ds:
        inside = ~read(ds.variables['study_mask']).astype(bool)
        n = int(inside.sum())
        for pol in POLLUTANTS:
            name = pol.replace('.', '')
            sums = np.zeros((3, n))
            counts = np.zeros(n, dtype=int)
            for day in range(len(ds.dimensions['time'])):
                layer = np.ma.filled(ds.variables[name + '_fraction'][day], np.nan)
                require(np.isnan(layer[:, ~inside]).all(), 'Outside fractions are not NaN')
                values = layer[:, inside]
                valid = np.isfinite(values).all(axis=0)
                require(np.array_equal(np.isfinite(values), np.broadcast_to(valid, values.shape)), 'Inconsistent source masks')
                require((values[:, valid] >= 0).all(), 'Negative fractions')
                require(np.allclose(values[2, valid], values[0, valid] + values[1, valid], rtol=3e-06, atol=1e-09), 'Daily fraction tier sum mismatch')
                sums[:, valid] += values[:, valid]
                counts += valid
            saved_counts = ds.variables[name + '_valid_days'][:][inside]
            require(np.array_equal(saved_counts, counts), 'Saved valid-day counts differ')
            average = np.full_like(sums, np.nan)
            np.divide(sums, counts[None, :], out=average, where=counts[None, :] > 0)
            saved = np.ma.filled(ds.variables[name + '_annual_mean_fraction'][:], np.nan)[:, inside]
            require(np.allclose(saved, average, rtol=3e-06, atol=1e-09, equal_nan=True), 'Saved annual daily-ratio means differ')
            ratio = np.ma.filled(ds.variables[name + '_annual_ratio'][:], np.nan)[:, inside]
            require(np.allclose(ratio[2], ratio[0] + ratio[1], rtol=3e-06, atol=1e-09, equal_nan=True), 'Annual matched-days ratio tier sum mismatch')

def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    log = OUTPUT / ('SOURCE_APPORTIONMENT_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '.jsonl')
    failed = 0
    with log.open('x') as handle:
        for year in YEARS:
            try:
                result = build_year(year)
            except Exception as exc:
                failed += 1
                result = dict(year=year, state='blocked', error=str(exc))
            handle.write(json.dumps(relative_metadata(result)) + '\n')
            handle.flush()
            print(year, result['state'], result.get('error', ''), flush=True)
    print('Share:', log, flush=True)
    require(failed == 0, 'Some years blocked; inspect log. Inputs and partial/final outputs retained.')
if __name__ == '__main__':
    main()
