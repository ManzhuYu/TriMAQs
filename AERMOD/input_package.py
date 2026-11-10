"""Build the final UOGD input package with consistent CSV tables and names."""
from pathlib import Path
from decimal import Decimal
from datetime import datetime
import calendar,csv,hashlib,json,math,re,shutil,sqlite3
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from pyproj import Transformer

CATEGORIES={
 'activity':('drilling_completion','Well_drilling_new_4','drilling_completion_data.csv'),
 'production':('production','Wellpad_production3','production_data.csv'),
 'compressors':('compressors','Compressors_4','compressor_data_with_elevation.csv'),
 'other_uogd':('other_tier2','ONG_related_4','other_tier2_data_with_elevation.csv'),
}
CASE_FIELDS=['source_category','source_id','year','pollutant','location_id','release_set_id','meteorological_grid_id','receptor_grid_ids','annual_emissions_tonnes']
SITE_FIELDS=['source_id','location_id','latitude','longitude','elevation_m','location_note']
RELEASE_FIELDS=['release_set_id','source_category','source_id','release_id','source_type','weight','latitude','longitude','easting_m','northing_m','elevation_m','height_m','temperature_K','velocity_m_s','diameter_m','x_length_m','y_length_m','angle_degrees','reference_year','reference_release_id','reference_pollutant','scc_codes','parameter_method']
MONTH_FIELDS=[calendar.month_abbr[m].lower()+'_tonnes' for m in range(1,13)]
DATE_FIELDS={}


def require(ok,message):
    if not ok:raise ValueError(message)


def ident(value):
    text=str(value).strip()
    try:return str(int(Decimal(text)))
    except Exception:return text


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def code_id(prefix,value):
    return prefix+hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()[:20]


def records(root):
    for database in sorted(root.glob('*/*/audit.sqlite3')):
        with sqlite3.connect(database.as_uri()+'?mode=ro',uri=True) as connection:
            for source,encoded in connection.execute('SELECT source,record FROM cases ORDER BY source'):
                yield database.parent/source,json.loads(encoded)
    for filename in sorted(root.glob('*/*/*/AERMOD_DONE.json')):
        db=filename.parent.parent/'audit.sqlite3'
        if db.is_file():
            with sqlite3.connect(db.as_uri()+'?mode=ro',uri=True) as c:
                if c.execute('SELECT 1 FROM cases WHERE source=?',(filename.parent.name,)).fetchone():continue
        yield filename.parent,json.loads(filename.read_text())


def read_releases(folder,record,transform):
    files=list(folder.glob('*.inp'));require(len(files)==1,'Expected one model input: '+str(folder))
    block=re.search(r'(?ms)^SO STARTING\s*\n.*?^SO FINISHED',files[0].read_text())
    require(block is not None,'Missing source section: '+str(folder))
    rows={};settings={}
    reference={r['model_source_id']:r for r in record.get('source_release_parameters',[])}
    for line in block[0].splitlines():
        p=line.split()
        if p and p[0]=='SO':p=p[1:]
        if p and p[0]=='LOCATION':
            x,y,z=map(float,p[3:6]);lon,lat=transform.transform(x,y)
            rows[p[1]]=dict(release_id=p[1],source_type=p[2],easting_m=x,northing_m=y,elevation_m=z,latitude=lat,longitude=lon)
        if p and p[0]=='SRCPARAM':settings[p[1]]=list(map(float,p[2:]))
    require(rows,'No model releases: '+str(folder))
    for sid,r in rows.items():
        v=settings[sid];require(len(v)>=5,'Incomplete release parameters')
        keys=['height_m','temperature_K','velocity_m_s','diameter_m'] if r['source_type']=='POINT' else ['height_m','x_length_m','y_length_m','angle_degrees'] if r['source_type']=='AREA' else None
        require(keys is not None,'Unsupported release type: '+r['source_type'])
        r.update(zip(keys,v[1:5]));ref=reference.get(sid,{})
        r['weight']=float(ref['weight']) if 'weight' in ref else 1. if len(rows)==1 else None
        require(r['weight'] is not None,'Missing multi-release allocation weights')
        r['reference_year']=ref.get('reference_year','')
        r['parameter_method']=ref.get('method',ref.get('parameter_method','Source parameters used in the AERMOD simulation'))
        r['reference_release_id']=ref.get('release_id','')
        r['reference_pollutant']=ref.get('reference_release_pollutant','')
        r['scc_codes']=';'.join(ref.get('SCCs',[]))
        require(all(math.isfinite(float(r[k])) for k in ['weight','latitude','longitude','elevation_m','height_m']),'Nonfinite release parameter')
        require(r['weight']>=0 and r['height_m']>=0,'Invalid weight/release height')
        if r['source_type']=='POINT':require(r['diameter_m']>0 and r['temperature_K']>0 and r['velocity_m_s']>=0,'Invalid point release')
        else:require(r['x_length_m']>0 and r['y_length_m']>0,'Invalid area release')
    require(abs(sum(r['weight'] for r in rows.values())-1)<1e-10,'Release allocation is not one')
    return list(rows.values())


def mass_values(row,workflow,year,expected_annual):
    """Standardize to daily or monthly totals, preserving the simulated mass."""
    row={str(k).strip().lower():v for k,v in row.items()}
    if workflow=='activity':
        days=366 if calendar.isleap(year) else 365
        if year not in DATE_FIELDS:DATE_FIELDS[year]=[(datetime(year,1,1)+pd.Timedelta(days=i)).strftime('%Y-%m-%d') for i in range(days)]
        fields=DATE_FIELDS[year]
        values=np.array([0. if pd.isna(row.get(c)) else float(row[c]) for c in fields],dtype=float)
        basis='daily_total'
    else:
        raw=[]
        for m in range(1,13):
            label=calendar.month_abbr[m].lower()
            if workflow=='production':
                matches=[k for k in row if re.fullmatch(label+r'(?:_(?:'+str(year)[-2:]+'|'+str(year)+r'))?_production',k)]
            else:matches=[label] if label in row else []
            require(len(matches)<=1,'Ambiguous month: '+label)
            v=row.get(matches[0]) if matches else None
            raw.append(0. if v is None or pd.isna(v) else float(v))
        raw=np.asarray(raw,dtype=float);lengths=np.array([calendar.monthrange(year,m)[1] for m in range(1,13)])
        monthly=float(raw.sum());daily=float((raw*lengths).sum())
        if math.isclose(monthly,expected_annual,rel_tol=1e-9,abs_tol=0.):values=raw;basis='monthly_total'
        elif math.isclose(daily,expected_annual,rel_tol=1e-9,abs_tol=0.):values=raw*lengths;basis='monthly_mean_daily_rate'
        else:raise ValueError('Inventory time basis does not reproduce the completed annual emissions')
        fields=MONTH_FIELDS
    require(np.isfinite(values).all() and (values>=0).all(),'Invalid inventory mass')
    require(math.isclose(float(values.sum()),expected_annual,rel_tol=1e-9,abs_tol=0.),'Annual emission mass differs from completed model')
    return fields,values,basis


def write_csv(path,fields,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)


def copy_support(stage,archive,source_root):
    """Copy only the final spatial inputs, without prior reports or alternatives."""
    mappings=[
      ('03_support/grids/grid_proj_1km_studyarea.csv','supporting/grids/study_grid_1km.csv'),
      ('03_support/grids/grid_proj_10km_studyarea_updated.csv','supporting/grids/meteorological_grid_10km.csv'),
      ('03_support/grids/grid_mapping_1km_10km.csv','supporting/grids/grid_links.csv'),
      ('03_support/current_mapping/receptor_to_grid_mapping_current.csv','supporting/receptor_mapping/receptor_mapping.csv'),
      ('03_support/current_mapping/reference_receptors_current.csv','supporting/receptor_mapping/reference_receptors.csv'),
      ('03_support/ml_grid/grid_spec_and_mask.npz','supporting/ml_grid/grid_spec_and_mask.npz')]
    for source_name,target_name in mappings:
        source=archive/source_name;require(source.is_file(),'Missing supporting input: '+str(source))
        target=stage/target_name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        require(sha(source)==sha(target),'Supporting copy differs')
    boundary=archive/'03_support/boundary'
    for suffix in ['.shp','.shx','.dbf','.prj','.cpg','.sbn','.sbx','.shp.xml']:
        source=boundary/('Counties_in_study_area'+suffix)
        if suffix in ['.shp','.shx','.dbf','.prj']:require(source.is_file(),'Missing boundary companion: '+str(source))
        if source.is_file():
            target=stage/'supporting/boundary'/('study_boundary'+suffix);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
            require(sha(source)==sha(target),'Boundary copy differs')
    files=sorted((archive/'03_support/current_receptors').glob('Grid_*/receptor_output_*_grid.DAT'))
    require(files,'No current receptor files')
    for source in files:
        target=stage/'supporting/receptors'/source.parent.name/source.name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        require(sha(source)==sha(target),'Receptor copy differs')
    # The old mapping file is not a second input: it provides only parent-grid links.
    import numpy as np
    grid=pd.read_csv(stage/'supporting/grids/study_grid_1km.csv')
    links=pd.read_csv(stage/'supporting/grids/grid_links.csv')
    require(len(files)==len(pd.read_csv(stage/'supporting/grids/meteorological_grid_10km.csv')),'Receptor grid count differs from the study grid')
    with np.load(stage/'supporting/ml_grid/grid_spec_and_mask.npz',allow_pickle=False) as spec:
        ids=grid.grid_ID.to_numpy(dtype=int);mask=spec['study_mask']
        require(set(ids)==set(np.flatnonzero(~mask.ravel())),'Study mask/grid IDs differ')
        require(np.allclose(spec['lat'].ravel()[ids],grid.lat,rtol=0,atol=5e-6) and np.allclose(spec['lon'].ravel()[ids],grid.lon,rtol=0,atol=5e-6),'Grid coordinates do not align')
    return len(files)


def make_readmes(stage):
    notes={
     'README.md':'''# UOGD model inputs

This folder contains the final emission inventories, source parameters and spatial inputs for the TriMAQs UOGD simulations, 2004–2023. Paths are relative to this folder. The accompanying code uses one dataset root and reads these files directly.

`emissions/` contains separate drilling/completion, production, compressor and other Tier 2 inventories. `parameters/` contains the final source locations, stack parameters and model case table. `supporting/` contains grids, the study boundary, receptors, spatial mappings and the ML grid mask. Each folder has a README. `data_dictionary.csv` defines the common fields, and `file_manifest.csv` lists file sizes and SHA256 checksums.

Tier 1 consists of drilling/completion and production. Tier 2 consists of compressors and other facilities related to refining, storage, gathering and transmission. The `other_tier2` group excludes the separately modeled compressors.

All shared inventory values are metric tonnes: daily totals for drilling/completion and monthly totals for the other categories. Allocation within a month is uniform by day. No additional component inventory should be added to production: leakage and heaters are already combined. NOX inventories are used for separate NOX and NO2 simulations, not added together. Missing inventory years reflect the modeled scope; production SO2 and other Tier 2 CH4 were not estimated.

The inventories contain the sources selected for the modeled study area. Source and release parameters retain the settings used in the simulations; they are not all measured values. Estimated geometry and reference-year assumptions are identified where recorded in the parameter table. Meteorological files, raw ACS data, upstream worksheets and earlier alternative inputs are not included here.
''',
     'parameters/README.md':'''# Source locations and model parameters

There is one final location table for each source category: `drilling_completion_data.csv`, `production_data.csv`, `compressor_data_with_elevation.csv`, and `other_tier2_data_with_elevation.csv`. Latitude and longitude are in degrees; elevation is metres above mean sea level. A location_id identifies each coordinate/elevation combination. If a source is represented at different locations in different cases, those entries remain explicit in the same table; the case table identifies which is used.

`stack_parameters.csv` contains point-source stack parameters (height, temperature, exit velocity and diameter), area-source dimensions and allocation weights. The source_type field identifies which fields apply. A release_set_id groups the releases that make up a modeled source. Weights sum to one for each set. Point rates are g/s; area rates are g/s/m2 after division by the rectangle area. Distances are metres, temperature is kelvin, exit velocity is m/s, and angle is degrees. Blank point/area-specific fields mean not applicable. Parameter_method and reference_year describe the basis where available; blank reference years should not be interpreted as measurements from the simulation year.

`model_cases.csv` links each source, year and pollutant to one location, release set, meteorological grid and receptor-grid list. A case is one source/pollutant/year simulation. The list of receptor grid IDs uses semicolons. These are ordinary model settings, read in the same way throughout the batch. Annual_emissions_tonnes is the integrated annual mass used to check the inventory and hourly allocation.

The shared code reads these CSVs directly. Location_note records source-assignment limitations where applicable. Estimated values remain estimates; the table names do not imply that every parameter was measured.
''',
     'supporting/README.md':'''# Spatial inputs

`grids/` contains the 1 km study grid, 10 km meteorological grid and links between them. `boundary/` contains the study-area shapefile and its companion files. `receptors/` contains the current AERMAP receptor definitions, grouped by 10 km grid. These files retain their model-readable names.

`receptor_mapping/receptor_mapping.csv` selects four receptors per study cell; a receptor may be selected for more than one cell. `reference_receptors.csv` lists the full coordinate reference used when summing source contributions. This is a four-receptor spatial estimate, not an exact area-weighted grid mean.

`ml_grid/grid_spec_and_mask.npz` gives the collaborator raster. x/y are column/row indices; lat/lon are geographic coordinates. True in study_mask means outside the study area. Study grid_ID is the zero-based flattened raster index. Concentrations outside the study area remain missing.
'''}
    notes.update({
     'supporting/grids/README.md':'# Study grids\n\nstudy_grid_1km.csv contains the study cells and their geographic coordinates. Grid_ID is the zero-based flattened index in the collaborator raster. meteorological_grid_10km.csv supplies the centers used for meteorology and receptor selection. grid_links.csv links each 1 km cell to a 10 km grid. Projected x/y values in the CSVs retain the source grid convention; use lat/lon and grid IDs for alignment with the NetCDF data. They are different from the raster-index x/y arrays in the NPZ file.\n',
     'supporting/boundary/README.md':'# Study boundary\n\nstudy_boundary.shp and its companion files define the study area. Keep the shp, shx, dbf and prj files together. The boundary is used to select modeled source origins; it is not the same as an individual source simulation footprint.\n',
     'supporting/receptors/README.md':'# Receptors\n\nThe Grid_ID folders contain the AERMAP receptor definitions used by AERMOD. File names are retained in the model-readable format receptor_output_ID_grid.DAT. Cartesian coordinates and terrain/height fields are in metres using the supplied UTM zone 18 convention. These are model input files rather than concentration outputs.\n',
     'supporting/receptor_mapping/README.md':'# Receptor sampling map\n\nreceptor_mapping.csv contains receptor easting/northing and the 1 km grid ID receiving that sample. Each study grid has four selected receptors; one receptor may serve more than one grid. reference_receptors.csv contains the complete coordinate list. Source contributions are summed at each unique coordinate before averaging the four samples per grid.\n',
     'supporting/ml_grid/README.md':'# Collaborator grid\n\ngrid_spec_and_mask.npz contains x, y, lat, lon, study_mask and crs. The raster has 815 rows and 842 columns. x/y are column/row indices, not metre coordinates. lat/lon are two-dimensional geographic grid-center coordinates. True in study_mask is outside the study area; False identifies the 222,897 study cells. The supplied CRS label is general; use the geographic coordinates for alignment.\n'
    })
    for category in ['drilling_completion','production','compressors','other_tier2']:
        prefix='compressor' if category=='compressors' else category
        timing='daily' if category=='drilling_completion' else 'monthly'
        notes[f'emissions/{category}/README.md']=f'''# {category.replace('_',' ').capitalize()} emission inventory

Files are grouped by pollutant and named `{prefix}_emissions_POLLUTANT_YEAR.csv`, with a lowercase pollutant in the filename (for example nox or pm25). Source_id is the source identifier and is read as text. Files contain {timing} total mass in metric tonnes. Daily fields are calendar dates; monthly fields are jan_tonnes through dec_tonnes. Source positions and stack parameters are in ../../parameters/.

Monthly totals are divided by the number of days in the month, then by 86400 seconds, to obtain g/s after multiplying tonnes by 1,000,000. Daily totals use the same conversion without division by days per month. Zeros indicate no emission in that period. The tables have complete numeric fields; missing required values should not be silently replaced during modeling. Annual totals are reconciled with the mass used in the completed simulations.

NOX mass feeds both NOX and NO2 model options. These are separate concentration outputs, not additive emission inventories. The {category.replace('_',' ')} category is modeled separately and combined with the appropriate tier during post-processing.
'''
    for name,text in notes.items():
        path=stage/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text,encoding='utf-8')
    fields=[
     ('source_id','Source identifier; read as text','1'),('source_category','drilling_completion, production, compressors or other_tier2','1'),
     ('year','Calendar year','year'),('pollutant','Modeled pollutant; NO2 uses the NOX inventory','1'),
     ('YYYY-MM-DD','Daily total emission mass','metric tonnes'),('jan_tonnes … dec_tonnes','Monthly total emission mass','metric tonnes'),
     ('latitude','Source or release latitude','degrees_north'),('longitude','Source or release longitude','degrees_east'),
     ('elevation_m','Ground elevation above mean sea level','m'),('location_id','Key to the location used for a case','1'),('location_note','Source-assignment note where applicable','text'),
     ('release_set_id','Key grouping releases and their allocation weights','1'),('release_id','AERMOD source identifier within a release set','1'),
     ('source_type','POINT or AREA','1'),('weight','Fraction of source emissions assigned to a release; sums to one within a set','1'),
     ('easting_m','Release easting in UTM zone 18N, as used by AERMOD','m'),('northing_m','Release northing in UTM zone 18N, as used by AERMOD','m'),
     ('height_m','Release height above local ground','m'),('temperature_K','Point exit-gas temperature','K'),('velocity_m_s','Point exit velocity','m/s'),
     ('diameter_m','Point stack diameter','m'),('x_length_m','Area rectangle x length','m'),('y_length_m','Area rectangle y length','m'),('angle_degrees','Area rectangle orientation','degree'),
     ('reference_year','Reference year for parameters where recorded','year'),('reference_release_id','Identifier of the reference inventory release where recorded','text'),
     ('reference_pollutant','Pollutant used for reference allocation where recorded','text'),('scc_codes','Semicolon-separated source classification codes where recorded','text'),
     ('parameter_method','Parameter source or estimation method where recorded','text'),
     ('meteorological_grid_id','10 km grid supplying AERMET meteorology','1'),('receptor_grid_ids','Semicolon-separated grid IDs supplying receptors','1'),
     ('annual_emissions_tonnes','Annual integrated source mass for the modeled case','metric tonnes')]
    write_csv(stage/'data_dictionary.csv',['field','description','units'],[dict(zip(['field','description','units'],r)) for r in fields])


def prepare_package(model_runs,archive,destination,boundary_file,compilation_summary,report_parent,include_support=True):
    model_runs,archive,destination=map(Path,[model_runs,archive,destination]);report_parent=Path(report_parent)
    report_parent.mkdir(parents=True,exist_ok=True)
    report_file=report_parent/('UOGD_INPUT_PACKAGE_SUMMARY_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')+'.json')
    stage=destination.parent/'working/input_packaging'
    require(not destination.exists(),'Final input folder already exists; it will not be overwritten')
    require(not stage.exists(),'An incomplete input build exists: '+str(stage))
    require(model_runs.is_dir(),'Completed model folder not found')
    boundary=gpd.read_file(boundary_file);require(boundary.crs is not None,'Boundary CRS missing')
    boundary=boundary.to_crs('EPSG:4326');domain=boundary.geometry.union_all() if hasattr(boundary.geometry,'union_all') else boundary.geometry.unary_union
    transform=Transformer.from_crs('EPSG:32618','EPSG:4326',always_xy=True)
    stage.mkdir(parents=True,exist_ok=False);sites={};release_sets={};selected={};seen=set();counts={};excluded=0;geometry_cache={};basis_counts={};mass_error=0.
    report=dict(state='blocked',destination=str(destination),errors=[])
    try:
        expected_rows=[];reviews={}
        aliases={v[1]:v[0] for v in CATEGORIES.values()};aliases['other_uogd']='other_tier2'
        if compilation_summary is not None:
            expected_rows=[json.loads(line) for line in Path(compilation_summary).read_text().splitlines() if line.strip()]
            require(all(r['state'] in {'completed','completed_with_review_flags','inventory_absent'} for r in expected_rows),'Compilation summary contains unfinished groups')
            for row in expected_rows:
                category=aliases.get(row['workflow'],row['workflow'])
                reviews.setdefault(category,set()).update(map(str,row.get('unresolved_source_assignment_review',[])))
        casefile=stage/'parameters/model_cases.csv';casefile.parent.mkdir(parents=True,exist_ok=True)
        with casefile.open('w',newline='',encoding='utf-8') as handle:
            writer=csv.DictWriter(handle,fieldnames=CASE_FIELDS);writer.writeheader()
            for workflow,(category,old_folder,sitefile) in CATEGORIES.items():
                folder=model_runs/old_folder if (model_runs/old_folder).is_dir() else model_runs/category
                require(folder.is_dir(),'Missing category model folder: '+str(folder))
                cases_read=0
                for casefolder,record in records(folder):
                    if record.get('state')!='completed':continue
                    require(record.get('model',{}).get('fatal_errors')==0,'Completed record lacks successful model status')
                    source=str(record['source_id']);year=int(record['year']);pol=record['pollutant']
                    require(2004<=year<=2023,'Case outside study period')
                    signature=(category,source,record['fingerprint'].get('script_sha256'),json.dumps(dict(releases=record.get('source_release_parameters',[]),stack=record.get('stack_parameters',{}),location=record['fingerprint'].get('coordinate_correction')),sort_keys=True))
                    if signature not in geometry_cache:geometry_cache[signature]=read_releases(casefolder,record,transform)
                    release=geometry_cache[signature]
                    within=[domain.covers(Point(r['longitude'],r['latitude'])) for r in release]
                    if not any(within):excluded+=1;continue
                    require(all(within),'Mixed inside/outside source components need review')
                    annual=float(record['annual_emissions_tonnes']);require(math.isfinite(annual) and annual>0,'Invalid annual case emissions')
                    origin=record['fingerprint'].get('coordinate_correction');stack=record.get('stack_parameters',{})
                    if origin:lat,lon=float(origin['latitude']),float(origin['longitude'])
                    elif 'easting' in stack:lon,lat=transform.transform(stack['easting'],stack['northing'])
                    else:lat,lon=release[0]['latitude'],release[0]['longitude']
                    location=dict(source_id=source,latitude=float(lat),longitude=float(lon),elevation_m=release[0]['elevation_m'],
                                  location_note='Source assignment requires review' if source in reviews.get(category,set()) else '')
                    locid=code_id('loc_',dict(category=category,**location));sites[(category,locid)]=dict(location,location_id=locid)
                    release_key=code_id('rel_',dict(category=category,source_id=source,releases=release))
                    release_sets[release_key]=(category,source,release)
                    grids=[]
                    for item in record['receptor_input_files']:
                        match=re.search(r'Grid_(\d+)',item['path']);require(match is not None,'Missing receptor grid ID');grids.append(match[1])
                    require(grids,'No receptor grids for case')
                    key=(category,year,pol,source);require(key not in seen,'Duplicate completed case');seen.add(key)
                    writer.writerow(dict(source_category=category,source_id=source,year=year,pollutant=pol,location_id=locid,release_set_id=release_key,
                        meteorological_grid_id=str(record['met_context']['grid_id']),receptor_grid_ids=';'.join(sorted(set(grids))),annual_emissions_tonnes=annual))
                    group=(category,year,pol);counts[group]=counts.get(group,0)+1
                    inventory=record['fingerprint']['inventory'];inputfile=Path(inventory['path'])
                    rawkey=(workflow,year,str(inputfile));data=selected.setdefault(rawkey,dict(expected=inventory,sources={}))
                    old=data['sources'].get(source)
                    require(old is None or math.isclose(old,annual,rel_tol=1e-10,abs_tol=0.),'NOX/NO2 annual mass differs for shared inventory')
                    data['sources'][source]=annual;cases_read+=1
                    if cases_read%10000==0:print(category,'cases read:',cases_read,flush=True)
                print(category,'selected cases:',cases_read,flush=True)
        require(seen,'No completed study cases selected')
        if compilation_summary is not None:
            for row in expected_rows:
                if row['state'] not in {'completed','completed_with_review_flags'}:continue
                category=aliases.get(row['workflow'],row['workflow']);key=(category,int(row['year']),row['pollutant'])
                require(counts.get(key,0)==int(row['included_sources']),'Case count differs from compiled results: '+str(key))
        if compilation_summary is not None:
            expected_keys={(aliases.get(row['workflow'],row['workflow']),int(row['year']),row['pollutant']) for row in expected_rows if row['state'] in {'completed','completed_with_review_flags'}}
            require(set(counts).issubset(expected_keys),'Completed cases occur outside compiled study groups')
        inventory_rows={}
        for (workflow,year,filename),selection in selected.items():
            category=CATEGORIES[workflow][0];original=Path(filename);stat=original.stat();expected=selection['expected']
            require(stat.st_size==expected['size'] and stat.st_mtime_ns==expected['mtime_ns'],'Working inventory changed since modeling: '+str(original))
            frame=pd.read_csv(original,dtype=str);lower={str(c).strip().lower():c for c in frame.columns}
            idcol=lower['well_id' if workflow=='activity' else 'wellpadid' if workflow=='production' else 'facility_id']
            ids=frame[idcol].map(ident)
            if workflow=='activity':ids=ids.str.zfill(10)
            require(not ids.duplicated().any(),'Duplicate inventory source identifiers')
            inventory_pol=None
            # The original filename identifies production pollutant; other folders identify it.
            if workflow=='production':
                match=re.fullmatch(r'Emission_(.+)_'+str(year)+r'_WellPad.csv',original.name);require(match is not None,'Unknown production inventory name');inventory_pol=match[1]
            else:inventory_pol=original.parent.name
            inventory_pol='NOX' if inventory_pol.upper() in {'NOX','NO2'} else inventory_pol.upper()
            token='pm25' if inventory_pol=='PM2.5' else inventory_pol.lower()
            prefix='compressor' if category=='compressors' else category
            target=Path('emissions')/category/inventory_pol/f'{prefix}_emissions_{token}_{year}.csv'
            rows=inventory_rows.setdefault(str(target),dict(fields=None,rows={}))
            require(set(selection['sources']).issubset(set(ids)),'Completed source missing from its inventory: '+str(original))
            for index,source in ids.items():
                if source not in selection['sources']:continue
                annual=selection['sources'][source];fields,values,basis=mass_values(frame.loc[index].to_dict(),workflow,year,annual)
                basis_counts[category+'|'+basis]=basis_counts.get(category+'|'+basis,0)+1
                mass_error=max(mass_error,abs(float(values.sum())-annual))
                entry=dict(source_id=source,**dict(zip(fields,values)))
                if source in rows['rows']:
                    require(np.allclose([rows['rows'][source][k] for k in fields],values,rtol=1e-10,atol=0.),'Conflicting inventories map to one final filename')
                else:rows['rows'][source]=entry
                rows['fields']=['source_id']+fields
        for name,table in inventory_rows.items():write_csv(stage/name,table['fields'],[table['rows'][key] for key in sorted(table['rows'])])
        for _,(category,_,name) in CATEGORIES.items():write_csv(stage/'parameters'/name,SITE_FIELDS,[row for (c,k),row in sorted(sites.items()) if c==category])
        release_rows=[]
        for key,(category,source,releases) in sorted(release_sets.items()):
            for item in releases:release_rows.append(dict(release_set_id=key,source_category=category,source_id=source,**item))
        write_csv(stage/'parameters/stack_parameters.csv',RELEASE_FIELDS,release_rows)
        receptor_files=copy_support(stage,archive,None) if include_support else 0
        make_readmes(stage)
        files=sorted(p for p in stage.rglob('*') if p.is_file())
        manifest=[dict(file=p.relative_to(stage).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in files]
        write_csv(stage/'file_manifest.csv',['file','bytes','sha256'],manifest)
        # No processing reports or previous alternatives are placed in the data package.
        require(not list(stage.rglob('*.json')),'Unexpected JSON input in final package')
        destination.parent.mkdir(parents=True,exist_ok=True);stage.rename(destination)
        report.update(state='completed',model_cases=len(seen),release_sets=len(release_sets),inventory_files=len(inventory_rows),location_rows=len(sites),
            source_counts={category:len({r['source_id'] for (c,k),r in sites.items() if c==category}) for category,_,_ in CATEGORIES.values()},
            case_counts=[dict(category=k[0],year=k[1],pollutant=k[2],cases=v) for k,v in sorted(counts.items())],
            outside_boundary_cases_excluded=excluded,source_geometry_files_read=len(geometry_cache),receptor_files=receptor_files,
            original_time_bases=basis_counts,maximum_annual_mass_difference_tonnes=mass_error,files=len(manifest)+1,
            total_bytes=sum(r['bytes'] for r in manifest),limits='Parameters and inventories reproduce the completed case settings. Underlying parameter/source-assignment uncertainties remain; a clean name is not new scientific validation.')
    except Exception as exc:
        report['errors'].append(str(exc));raise
    finally:
        report_file.write_text(json.dumps(report,indent=2));print('Build summary (outside sharing folder):',report_file,flush=True)
    print('Final input package:',destination,flush=True)
    return report
