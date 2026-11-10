"""Prepare terrain and receptor inputs with AERMAP."""
from pathlib import Path as _CodePath
import sys as _code_sys
_code_sys.path.insert(0, str(_CodePath(__file__).resolve().parents[2]))
from trimaqs_config import path as _configured_path, setting as _setting
import os
import re
import shutil
import subprocess
import pandas as pd
import numpy as np
import logging
from shapely.geometry import Point, box
from pyproj import CRS, Transformer
import rasterio
import geopandas as gpd
from concurrent.futures import ThreadPoolExecutor
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=[logging.StreamHandler()])
output_base_dir_template = _configured_path('run_aermap.output_folder')
aermap_executable = _configured_path('run_aermap.executable')
dem_dir = _configured_path('run_aermap.dem_folder')
temp_dem_dir = _configured_path('run_aermap.temporary_dem_folder')
os.makedirs(temp_dem_dir, exist_ok=True)

def sanitize_filename(name):
    """Remove or replace illegal characters from filenames."""
    invalid_chars = '[<>:"/\\\\|?*]'
    return re.sub(invalid_chars, '_', name)

def convert_to_utm(lat, lon, utm_zone):
    """Convert geographic coordinates (lat/lon) to UTM coordinates."""
    utm_crs = CRS.from_dict({'proj': 'utm', 'zone': utm_zone, 'datum': 'NAD83', 'units': 'm', 'no_defs': True})
    transformer = Transformer.from_crs('epsg:4269', utm_crs, always_xy=True)
    return transformer.transform(lon, lat) + (utm_zone,)

def define_domain_bounds_utm(easting, northing, half_side_m=32000, uniform_margin_m=5000):
    """Define modeling domain bounds in UTM coordinates."""
    return {'min_easting': easting - half_side_m + uniform_margin_m, 'max_easting': easting + half_side_m - uniform_margin_m, 'min_northing': northing - half_side_m + uniform_margin_m, 'max_northing': northing + half_side_m - uniform_margin_m}

def generate_receptors(center_x, center_y, grid_id):
    """Use the supplied receptor coordinates without defining a different grid."""
    from pathlib import Path
    source = Path(_configured_path('run_aermap.coordinate_templates')) / f'Grid_{grid_id}' / f'receptor_output_{grid_id}_grid.DAT'
    output = Path(output_base_dir_template) / f'Grid_{grid_id}'
    if output.resolve() == source.parent.resolve():
        raise ValueError('AERMAP output must differ from coordinate templates')
    receptors = []
    for line in source.read_text().splitlines():
        parts = line.split()
        if parts and parts[0] == 'RE':
            parts = parts[1:]
        if parts and parts[0] in {'EVALCART', 'DISCCART'}:
            receptors.append((float(parts[1]), float(parts[2]), 0.0, 2, str(grid_id)[:8], f'{grid_id}_{len(receptors) + 1}'))
    if not receptors:
        raise ValueError('No supplied receptor coordinates: ' + str(source))
    return receptors

def find_closest_dems(lat, lon, dem_files, dem_crs, source_crs):
    """Find the 4 nearest DEM files covering the target location."""
    dem_geometries = []
    for dem_file in dem_files:
        with rasterio.open(dem_file) as src:
            dem_geometries.append(box(*src.bounds))
    dems_gdf = gpd.GeoDataFrame({'geometry': dem_geometries}, crs=CRS.from_epsg(dem_crs))
    projected_crs = CRS.from_dict({'proj': 'utm', 'zone': 18, 'datum': 'NAD83', 'units': 'm', 'no_defs': True})
    dems_projected = dems_gdf.to_crs(projected_crs)
    source_x, source_y = Transformer.from_crs(source_crs, projected_crs, always_xy=True).transform(lon, lat)
    source_point_projected = Point(source_x, source_y)
    dems_projected['centroid'] = dems_projected.geometry.centroid
    dems_projected['distance'] = dems_projected['centroid'].distance(source_point_projected)
    closest_dems = dems_projected.sort_values('distance').head(4)
    domain_bounds_utm = define_domain_bounds_utm(source_x, source_y)
    corners = [(domain_bounds_utm['min_easting'], domain_bounds_utm['min_northing']), (domain_bounds_utm['min_easting'], domain_bounds_utm['max_northing']), (domain_bounds_utm['max_easting'], domain_bounds_utm['min_northing']), (domain_bounds_utm['max_easting'], domain_bounds_utm['max_northing'])]
    for corner in corners:
        if not any((dem.contains(Point(corner)) for dem in closest_dems.geometry)):
            raise ValueError('DEM coverage insufficient.')
    return [dem_files[idx] for idx in closest_dems.index]

def generate_aermap_inp(facility_id, unit_id, temp_dem_files, output_dir, easting, northing, utm_zone, domain_bounds_utm, receptor_dat_file):
    """Generate AERMAP input file (.inp) for a given grid cell."""
    inp_file_path = os.path.join(output_dir, f'aermap_{sanitize_filename(facility_id)}_{sanitize_filename(unit_id)}.inp')
    try:
        with open(inp_file_path, 'w') as f:
            f.write('CO STARTING\n')
            f.write(f'   TITLEONE  AERMAP Run for {facility_id}_{unit_id}\n\n')
            f.write('   DATATYPE  NED\n')
            f.write('   FLAGPOLE 2\n')
            for dem_file in temp_dem_files:
                f.write(f'   DATAFILE  "{dem_file}"\n')
            f.write(f"   DOMAINXY  {domain_bounds_utm['min_easting']:.1f} {domain_bounds_utm['min_northing']:.1f} {utm_zone} {domain_bounds_utm['max_easting']:.1f} {domain_bounds_utm['max_northing']:.1f} {utm_zone}\n")
            f.write(f'   ANCHORXY  {easting:.1f} {northing:.1f}  {easting:.1f} {northing:.1f} {utm_zone} 2\n\n')
            f.write('   DEBUGOPT  ALL\n')
            f.write('   RUNORNOT  RUN\n')
            f.write('CO FINISHED\n\n')
            f.write('SO STARTING\n')
            f.write(f'   LOCATION {facility_id}_{unit_id} POINT {easting:.1f} {northing:.1f}\n')
            f.write('SO FINISHED\n\n')
            f.write('RE STARTING\n')
            f.write(f'   INCLUDED  "{receptor_dat_file}"\n')
            f.write('RE FINISHED\n\n')
            f.write('OU STARTING\n')
            f.write(f'   SOURCLOC  "source_output_{facility_id}_{unit_id}.DAT"\n')
            f.write(f'   RECEPTOR  "receptor_output_{facility_id}_{unit_id}.DAT"\n')
            f.write('OU FINISHED\n')
        logging.info(f'AERMAP input file generated: {inp_file_path}')
    except Exception as e:
        logging.error(f'Failed to write AERMAP input file: {e}')
    return inp_file_path

def run_aermap_in_output_dir(aermap_exec_path, inp_file_path, output_dir):
    """Run AERMAP in the specified output directory."""
    try:
        result = subprocess.run([aermap_exec_path, inp_file_path], cwd=output_dir, capture_output=True, text=True, check=True)
        logging.info(f'AERMAP run completed in: {output_dir}')
    except subprocess.CalledProcessError as e:
        logging.error('AERMAP run failed.')
        logging.error(e.stdout)
        logging.error(e.stderr)
        raise

def get_dem_crs(dem_files):
    """Ensure all DEMs share a single CRS (EPSG code)."""
    crs_set = {rasterio.open(f).crs.to_epsg() for f in dem_files if rasterio.open(f).crs}
    return crs_set.pop() if len(crs_set) == 1 else None

def process_grid(row, dem_files, dem_crs, source_crs, output_base_dir, utm_zone):
    """Process a single grid: generate receptors, prepare DEMs, run AERMAP."""
    grid_id, lat, lon = (str(int(row['grid_ID'])), row['lat'], row['lon'])
    logging.info(f'Processing Grid ID: {grid_id}')
    center_x, center_y, _ = convert_to_utm(lat, lon, utm_zone)
    output_dir = os.path.join(output_base_dir, f'Grid_{grid_id}')
    os.makedirs(output_dir, exist_ok=True)
    receptor_file = os.path.join(output_dir, f'receptors_{grid_id}_18.dat')
    with open(receptor_file, 'w') as f:
        for rec in generate_receptors(center_x, center_y, grid_id):
            f.write(f'RE EVALCART {rec[0]:.8f} {rec[1]:.8f} {rec[2]:.1f} {rec[3]} {rec[4]} {rec[5]}\n')
    try:
        closest_dems = find_closest_dems(lat, lon, dem_files, dem_crs, source_crs)
    except Exception as e:
        raise ValueError(f'DEM selection failed for Grid {grid_id}: {e}')
    temp_dir = os.path.join(temp_dem_dir, f'{grid_id}')
    os.makedirs(temp_dir, exist_ok=True)
    temp_dem_files = [shutil.copy(dem, temp_dir) for dem in closest_dems]
    inp_file = generate_aermap_inp(grid_id, 'grid', temp_dem_files, output_dir, center_x, center_y, utm_zone, define_domain_bounds_utm(center_x, center_y), receptor_file)
    run_aermap_in_output_dir(aermap_executable, inp_file, output_dir)
    shutil.rmtree(temp_dir, ignore_errors=True)

def main():
    grid_csv = _configured_path('run_aermap.grid_file')
    os.makedirs(output_base_dir_template, exist_ok=True)
    try:
        grids_df = pd.read_csv(grid_csv)
        grids_df = grids_df.rename(columns={c: 'grid_ID' if c.lower() == 'grid_id' else 'lat' if c.lower() == 'center_lat' else 'lon' if c.lower() == 'center_lon' else c for c in grids_df.columns})
    except Exception as e:
        raise ValueError(f'Failed to read grid CSV: {e}')
    dem_files = [os.path.join(dem_dir, f) for f in os.listdir(dem_dir) if f.endswith('.tif')]
    dem_crs = get_dem_crs(dem_files)
    if not dem_crs or dem_crs != 4269:
        raise ValueError('DEM files must share EPSG:4269')
    logging.info(f'Processing {len(grids_df)} grids in parallel.')
    source_crs = CRS.from_epsg(4269)
    with ThreadPoolExecutor(max_workers=7) as executor:
        futures = [executor.submit(process_grid, row, dem_files, dem_crs, source_crs, output_base_dir_template, 18) for _, row in grids_df.iterrows()]
        for future in futures:
            future.result()
if __name__ == '__main__':
    main()
