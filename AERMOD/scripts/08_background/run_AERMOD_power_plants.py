"""Prepare and run Appendix C power-plant/other-industry cases in scratch.

Uses the current shared regional inventories, current AERMET output 2 and
AERMAP output 2. Imports output validation/storage helpers from the unchanged
run_AERMOD_full_calculation.py beside this file. Requires numpy, geopandas,
pyproj and rasterio in the HPC environment.
"""
from pathlib import Path
import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[_v]='1'
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from contextlib import ExitStack
from datetime import datetime,timedelta
from functools import lru_cache
import argparse,calendar,csv,json,math,re,sqlite3,subprocess,traceback,zlib
import numpy as np
import run_AERMOD_full_calculation as core

BASE=Path(os.environ.get('TRIMAQS_DATA_ROOT',str(Path(__file__).resolve().parents[3]/'TriMAQs_data')))
VERSION='2026-11-07-appendix-background-v2'
PREDECESSOR_SHA='3e904161a3af7566657f9bf42735ae724da81c6b5b5a94e1b5076c5ca5415680'
CFG={};SITES={};STACKS={};MEDIANS={};GRID=[];DOMAIN=None;PROJECT=None
FOLDERS={'power_plants':'Power_plants','other_industry':'NEI_industrial'}
POLS=['NO2','PM2.5','SO2','VOC']
BLANK={'','nan','none','na','n/a','null'}


def norm(text):return re.sub(r'[^a-z0-9]+','_',str(text).lower()).strip('_')
def ident(text):return re.sub(r'^(\d+)\.0$',r'\1',str(text).strip())


def rows(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):yield {norm(k):v for k,v in row.items()}


def number(value):
    if str(value).strip().lower() in BLANK:return None
    result=float(value);core.require(math.isfinite(result),'Nonfinite input')
    return result


def alias(row,names,required=True):
    values=[row[k] for k in names if k in row and str(row[k]).strip().lower() not in BLANK]
    core.require(not values or len(set(str(v).strip() for v in values))==1,'Conflicting input aliases: '+str(names))
    core.require(values or not required,'Missing input column/value: '+str(names))
    return values[0] if values else ''


@lru_cache(maxsize=4)
def inventory(path):
    path=Path(path);before=core.file_info(path);result={}
    for row in rows(path):
        facility=ident(alias(row,['facilityid','facility_id','source_id','eis_facility_id']))
        unit=ident(alias(row,['unitid','unit_id'],False)) if CFG['category']=='power_plants' else ''
        key=(facility,unit);core.require(key not in result,'Duplicate inventory facility/unit: '+str(key))
        result[key]=row
    core.require(core.file_info(path)==before,'Inventory changed while reading')
    return result,dict(before,sha256=core.sha(path))


def inventory_path(year,pol):
    root=Path(CFG['data_root'])/'inputs/emissions'
    if CFG['category']=='power_plants':
        token={'NO2':'nox','SO2':'so2','PM2.5':'PM25','VOC':'VOC'}[pol]
        return root/'power_plants'/f'Daily {year}'/f'formatted-daily-{year}-{token}.csv'
    return root/'other_industry'/f"{'NOX' if pol=='NO2' else pol}_interpolated_emission.csv"


def emissions(case):
    year=int(case['year']);n=366 if calendar.isleap(year) else 365
    days=[datetime(year,1,1)+timedelta(days=i) for i in range(n)]
    table,info=inventory(str(inventory_path(year,case['pollutant'])))
    row=table[(case['facility_id'],case['unit_id'])]
    if CFG['category']=='power_plants':
        expected={d.strftime('%Y_%m_%d') for d in days}
        actual={k for k in row if re.fullmatch(r'\d{4}_\d{2}_\d{2}',k)}
        core.require(actual==expected,'Power-plant daily columns do not cover the exact calendar')
        raw=[number(row[d.strftime('%Y_%m_%d')]) for d in days]
        core.require(all(v is not None and v>=0 for v in raw),'Missing/negative daily power-plant emissions')
        grams={'short_ton':907184.74,'metric_tonne':1e6}[CFG['powerplant_mass_unit']]
        mass=[v*grams/1e6 for v in raw]
    else:
        total=number(row[str(year)+'_emission']);core.require(total is not None and total>=0,'Missing/negative industrial annual mass')
        mass=[total/n]*n
    return days,mass,info,row


def numeric_unit(value):
    text=ident(value)
    return str(int(text)) if re.fullmatch(r'\d+',text) else text


def table_row(table,facility,unit,power):
    """Use an exact unit, or a unique numeric-padding alias at that facility."""
    key=(facility,unit)
    if key in table:return table[key]
    if not power:return {}
    matches=[r for (f,u),r in table.items() if f==facility and numeric_unit(u)==numeric_unit(unit)]
    core.require(len(matches)<=1,'Ambiguous numeric unit alias at facility '+facility+': '+unit)
    return matches[0] if matches else {}


def read_site_table(path,power,kind='stack'):
    table={}
    for r in rows(path):
        facility=ident(alias(r,['facilityid','facility_id','source_id']))
        unit=ident(alias(r,['unitid','unit_id'])) if power else ''
        key=(facility,unit)
        if key not in table:table[key]=r;continue
        old=table[key]
        if old==r:continue
        if kind=='elevation':
            candidates=old.get('_location_candidates',[old])
            if r not in candidates:candidates.append(r)
            table[key]={'_location_candidates':candidates}
        else:
            # SCC, names or other descriptive columns can differ without a
            # change to the legacy physical POINT representation. Distinct
            # release IDs or physical parameters still require review.
            def physical_signature(row):
                return (stack_raw(row,power),
                    alias(row,['emission_release_point_type','erptype'],False),
                    alias(row,['emission_release_point_id','rel_point_id'],False))
            try:agree='_stack_conflicts' not in old and physical_signature(old)==physical_signature(r)
            except (ValueError,KeyError):agree=False
            if not agree:
                candidates=old.get('_stack_conflicts',[old])
                if r not in candidates:candidates.append(r)
                table[key]={'_stack_conflicts':candidates}
    return table


def stack_raw(row,power):
    def get(names):return number(alias(row,names,False))
    height=get(['stack_height','stack_height_feet'])
    temp=get(['exit_gas_temperature','exit_temperature_100_fahrenheit'])
    velocity=get(['exit_gas_velocity','exit_velocity_100_feet_per_second'])
    diameter=get(['stack_diameter','stack_diameter_feet'])
    if power and diameter is None:
        area=get(['area_at_top_square_feet'])
        if area is not None and area>0:diameter=math.sqrt(4*area/math.pi)
    return dict(height_m=None if height is None else height*.3048,
        temperature_K=None if temp is None else (temp-32)*5/9+273.15,
        velocity_m_s=None if velocity is None else velocity*.3048,
        diameter_m=None if diameter is None else diameter*.3048)


def valid(key,value):
    if value is None or not math.isfinite(value):return False
    if key=='height_m':return value>=0
    if key=='temperature_K':return value>0
    if key=='velocity_m_s':return value>=0
    return value>0


def stack_for(case):
    power=CFG['category']=='power_plants'
    row=table_row(STACKS,case['facility_id'],case['unit_id'],power)
    if power:core.require(row,'Power-plant facility/unit stack match missing; no facility-only substitution')
    core.require('_stack_conflicts' not in row,'Conflicting stack/release rows for this facility/unit; release-level review required')
    release_type=alias(row,['emission_release_point_type','erptype'],False)
    if release_type:
        core.require(ident(release_type) in {'2','02'},'Nonvertical/unknown NEI release type needs a source-specific representation: '+release_type)
    raw=stack_raw(row,CFG['category']=='power_plants');filled=[]
    for name,value in raw.items():
        if value is None:
            core.require(name in MEDIANS,'No valid category median for '+name)
            raw[name]=MEDIANS[name];filled.append(name)
        else:core.require(valid(name,value),'Invalid existing stack value: '+name)
    # Tiny nonmissing values can represent fugitive sources; do not median-fill them.
    core.require(raw['diameter_m']>=.01 and raw['velocity_m_s']>=.001,
        'Tiny diameter/velocity requires release-type review; not automatically replaced with stack medians')
    return raw,filled,'Recorded vertical stack' if release_type else 'Historical POINT representation; release type not recorded'


def initialize(config):
    global CFG,SITES,STACKS,MEDIANS,GRID,DOMAIN,PROJECT
    CFG=config;core.CONFIG=config
    power=config['category']=='power_plants'
    SITES=read_site_table(config['elevation_file'],power,kind='elevation')
    STACKS=read_site_table(config['stack_file'],power)
    MEDIANS={}
    for key in ['height_m','temperature_K','velocity_m_s','diameter_m']:
        good=[]
        for row in STACKS.values():
            typ=alias(row,['emission_release_point_type','erptype'],False)
            if typ and ident(typ) not in {'2','02'}:continue
            try:value=stack_raw(row,power)[key]
            except (ValueError,KeyError):continue
            if valid(key,value) and (key!='diameter_m' or value>=.01) and (key!='velocity_m_s' or value>=.001):good.append(value)
        if good:MEDIANS[key]=float(np.median(good))
    import geopandas as gpd
    from pyproj import Transformer
    boundary=gpd.read_file(config['boundary'])
    core.require(boundary.crs is not None,'Study boundary has no CRS')
    boundary=boundary.to_crs(4326)
    DOMAIN=boundary.geometry.union_all() if hasattr(boundary.geometry,'union_all') else boundary.geometry.unary_union
    PROJECT=Transformer.from_crs(4326,32618,always_xy=True)
    GRID=list(rows(config['grid_file']))
    core.require(GRID,'Meteorological grid table empty')


def distances(lat,lon):
    a=np.radians([float(r['center_lat']) for r in GRID]);b=np.radians([float(r['center_lon']) for r in GRID])
    la=math.radians(lat);lo=math.radians(lon)
    hav=np.sin((a-la)/2)**2+np.cos(a)*math.cos(la)*np.sin((b-lo)/2)**2
    return 12742*np.arcsin(np.sqrt(np.clip(hav,0,1)))


@lru_cache(maxsize=20000)
def elevation(lat,lon):
    import rasterio
    from pyproj import Transformer
    matches=[]
    for path in sorted(Path(CFG['dem_root']).rglob('*.tif')):
        with rasterio.open(path) as dem:
            core.require(dem.crs is not None,'DEM has no CRS: '+str(path))
            x,y=Transformer.from_crs(4326,dem.crs,always_xy=True).transform(lon,lat)
            if not (dem.bounds.left<=x<dem.bounds.right and dem.bounds.bottom<y<=dem.bounds.top):continue
            core.require(not dem.units or dem.units[0] in {None,'','m','meter','metre','meters','metres'},'DEM vertical units are not metres')
            value=next(dem.sample([(x,y)],masked=True))[0]
            if not np.ma.is_masked(value) and math.isfinite(float(value)):matches.append(float(value))
    core.require(matches and max(matches)-min(matches)<=2,'Missing or conflicting DEM elevation at source coordinates')
    return float(np.median(matches))


def coordinates(case,row):
    site=table_row(SITES,case['facility_id'],case['unit_id'],CFG['category']=='power_plants')
    candidates=site.get('_location_candidates',[site]) if site else []
    # Inventory coordinates define the emission source. Elevation is reused only
    # if its table coordinates describe the same location (within ~10 m).
    lat=number(alias(row,['latitude','lat'],False));lon=number(alias(row,['longitude','lon'],False))
    if lat is None or lon is None:
        positions=set()
        for candidate in candidates:
            a=number(alias(candidate,['latitude','lat'],False));b=number(alias(candidate,['longitude','lon'],False))
            if a is not None and b is not None and (lat is None or abs(a-lat)<=.0001) and (lon is None or abs(b-lon)<=.0001):positions.add((a,b))
        core.require(len(positions)==1,'Source coordinates missing/ambiguous across elevation-table rows; need inventory coordinates')
        lat,lon=next(iter(positions))
    core.require(lat is not None and lon is not None and -90<=lat<=90 and -180<=lon<=180,'Invalid source coordinates')
    from shapely.geometry import Point
    if not DOMAIN.covers(Point(lon,lat)):return None
    z=None;heights=[]
    for candidate in candidates:
        slat=number(alias(candidate,['latitude','lat'],False));slon=number(alias(candidate,['longitude','lon'],False))
        if slat is not None and slon is not None and abs(slat-lat)<=.0001 and abs(slon-lon)<=.0001:
            value=number(alias(candidate,['elevation','elevation_m'],False))
            if value is not None:heights.append(value)
    if heights and max(heights)-min(heights)<=.01:z=float(np.median(heights))
    if z is None:z=elevation(lat,lon)
    x,y=PROJECT.transform(lon,lat)
    return dict(latitude=lat,longitude=lon,elevation_m=z,easting_m=x,northing_m=y)


@lru_cache(maxsize=12)
def met_context(year,gid):
    folder=Path(CFG['met_root'])/str(year)/f'Grid_{gid}'
    manifest=json.loads((folder/'manifest.json').read_text());audit=json.loads((folder/'output_audit.json').read_text())
    core.require(ident(manifest['grid_id'])==gid and int(manifest['year'])==year,'Meteorological manifest grid/year differs')
    core.require(audit['hours']==(8784 if calendar.isleap(year) else 8760),'Meteorological annual clock incomplete')
    bad=int(audit['noncalm_invalid_hours']);exceptions=[]
    if bad:
        core.require(bad==1 and len(audit['invalid_times'])==1,'Unreviewed meteorological exceptions')
        wanted=set(audit['invalid_times'])
        for line in (folder/'output_stage2.SFC').read_text().splitlines()[1:]:
            p=line.split()
            if not p:continue
            yy=int(p[0]);yy=yy+2000 if yy<100 else yy
            time=datetime(yy,int(p[1]),int(p[2]))+timedelta(hours=int(p[4]))
            if str(time) in wanted:
                core.require(float(p[5])==0 and float(p[7])==-9 and float(p[9])==-999 and float(p[11])==-8888,'Meteorological exception is not the reviewed neutral pattern')
                exceptions.append(dict(time_EST=str(time),SFC_record=line))
        core.require(len(exceptions)==bad,'Reviewed exception record missing')
    else:core.require((folder/'DONE.json').is_file() and manifest.get('state')=='completed_pending_scientific_review','Meteorology has no completion receipt')
    for name in ['AERMET_REPORT.RPT','stage2.rpt']:
        text=(folder/name).read_text(errors='replace').upper()
        core.require('AERMET FINISHED SUCCESSFULLY' in text and 'UN-SUCCESSFULLY' not in text,'AERMET stage not successful')
    files={n:core.file_info(folder/n) for n in ['output_stage2.SFC','output_stage2.PFL']}
    core.require(all(v['size']>0 for v in files.values()),'Empty met file')
    return dict(grid_id=gid,year=year,elevation_m=float(manifest['elevation_m']),exceptions=exceptions,files=files)


def receptor_network(lat,lon,radius,work):
    points={};files=[]
    for row,dist in zip(GRID,distances(lat,lon)):
        if dist>radius:continue
        gid=ident(row['grid_id']);path=Path(CFG['aermap_root'])/f'Grid_{gid}'/f'receptor_output_{gid}_grid.DAT'
        files.append(core.file_info(path))
        for line in path.read_text().splitlines():
            p=core.fields(line)
            if not p or p[0] not in {'EVALCART','DISCCART'}:continue
            xy=(round(float(p[1]),2),round(float(p[2]),2))
            core.require(xy not in points or points[xy][1:6]==p[1:6],'Conflicting receptor terrain/flagpole')
            points[xy]=p
    core.require(points,'No receptor-grid centres within the selected radius')
    path=work/'receptors.DAT'
    path.write_text('RE ELEVUNIT METERS\n'+'\n'.join('RE '+' '.join(p) for p in points.values())+'\n',encoding='ascii')
    return set(points),files


def control(case,release,context):
    pol=case['pollutant'];year=case['year'];sid=release['release_id']
    options={'NO2':'CONC ALPHA WETDPLT TTRM2 ARM2','PM2.5':'CONC ALPHA WETDPLT DRYDPLT','SO2':'CONC ALPHA WETDPLT','VOC':'CONC'}[pol]
    lines=['CO STARTING','   TITLEONE Appendix C '+case['source_id'],'   MODELOPT '+options]
    if pol=='NO2':lines+=['   OZONEVAL 40.0 UG/M3','   NO2STACK 0.10']
    lines+=['   AVERTIME 24','   POLLUTID '+pol,'   RUNORNOT RUN','   ERRORFIL "model.err"','CO FINISHED',
        'SO STARTING','   ELEVUNIT METERS',f"   LOCATION {sid} POINT {release['easting_m']:.8f} {release['northing_m']:.8f} {release['elevation_m']:.8f}",
        f"   SRCPARAM {sid} 0 {release['height_m']:.8f} {release['temperature_K']:.8f} {release['velocity_m_s']:.8f} {release['diameter_m']:.8f}",
        '   HOUREMIS "hourly_emissions.txt" '+sid]
    if pol=='PM2.5':lines += [f'   PARTDIAM {sid} 0.1',f'   PARTDENS {sid} 1.0',f'   MASSFRAX {sid} 1.0']
    if pol in {'NO2','SO2'}:lines += [f"   GASDEPOS {sid} "+('0.1361 1.4E-5 12000 8444' if pol=='NO2' else '0.1122 1.83E-5 732 72')]
    lines+=['   SRCGROUP ALL','SO FINISHED','RE STARTING','   INCLUDED "receptors.DAT"','RE FINISHED','ME STARTING',
        '   SURFFILE "'+context['files']['output_stage2.SFC']['path']+'"',
        '   PROFFILE "'+context['files']['output_stage2.PFL']['path']+'"',f'   SURFDATA 99999 {year} STUDY_AREA',
        f'   UAIRDATA 99999 {year} STUDY_AREA',f"   PROFBASE {context['elevation_m']:.8f}",
        f'   STARTEND {year} 01 01 {year} 12 31','ME FINISHED','OU STARTING','   SUMMFILE "model.sum"','   DAYTABLE 24','OU FINISHED']
    return '\n'.join(lines)+'\n'


def radius_km(annual_metric_tonnes):
    short_tons=annual_metric_tonnes/.90718474
    # Round only conversion noise at exact thresholds, not physical values.
    if math.isclose(short_tons,3600,rel_tol=1e-12):short_tons=3600.
    if math.isclose(short_tons,1000,rel_tol=1e-12):short_tons=1000.
    return 30 if short_tons>=3600 else 20 if short_tons>=1000 else 15


def radius_fingerprints(payload):
    current=core.signature(payload)
    previous=dict(payload,config=dict(payload['config'],script_sha256=PREDECESSOR_SHA))
    return current,core.signature(previous)


def commit_result(group,sid,record,work):
    with core.database(group) as conn:
        conn.execute('CREATE TABLE IF NOT EXISTS radius_history (source TEXT,fingerprint TEXT,record TEXT,PRIMARY KEY(source,fingerprint))')
        conn.execute('CREATE TABLE IF NOT EXISTS radius_history_diagnostics (source TEXT,fingerprint TEXT,name TEXT,content BLOB,PRIMARY KEY(source,fingerprint,name))')
        old=conn.execute('SELECT record FROM cases WHERE source=?',(sid,)).fetchone()
        if old and json.loads(old[0])['fingerprint']!=record['fingerprint']:
            prior=json.loads(old[0]);core.require(prior['fingerprint']==record.get('previous_radius_fingerprint'),'Unexpected previous completion while committing')
            prior['previous_canonical_csv']=prior['output_files'][0]['path']
            prior['output_files'][0]['path']=record['previous_radius_csv']
            conn.execute('INSERT OR IGNORE INTO radius_history VALUES (?,?,?)',(sid,prior['fingerprint'],json.dumps(prior)))
            for name,content in conn.execute('SELECT name,content FROM diagnostics WHERE source=?',(sid,)).fetchall():
                conn.execute('INSERT OR IGNORE INTO radius_history_diagnostics VALUES (?,?,?,?)',(sid,prior['fingerprint'],name,content))
        conn.execute('INSERT OR REPLACE INTO cases VALUES (?,?)',(sid,json.dumps(record)))
        conn.execute('DELETE FROM diagnostics WHERE source=?',(sid,))
        for name in ['aermod.inp','model.err','model.sum','run.log']:
            if (work/name).is_file():conn.execute('INSERT INTO diagnostics VALUES (?,?,?)',(sid,name,zlib.compress((work/name).read_bytes())))


def recover_radius_update(group,sid):
    work=group/(sid+'.running');receipt=work/'receipt.json';target=group/f'aermod_{sid}.csv.gz'
    if not receipt.is_file():return
    record=json.loads(receipt.read_text());info=record['output_files'][0]
    core.require(record['source_id']==sid and info['path']==str(target),'Unexpected interrupted receipt')
    if target.is_file() and core.sha(target)==info['sha256']:
        commit_result(group,sid,record,work);return
    candidate=work/target.name
    core.require(candidate.is_file() and core.sha(candidate)==info['sha256'],'Interrupted revised result missing/changed')
    if target.exists():
        prior=record.get('previous_radius_output')
        core.require(prior and core.sha(target)==prior['sha256'],'Unexpected previous CSV during recovery')
        backup=Path(record['previous_radius_csv']);backup.parent.mkdir(parents=True,exist_ok=True)
        core.require(not backup.exists(),'Backup already exists with a previous CSV still at target')
        target.rename(backup)
    candidate.rename(target);commit_result(group,sid,record,work)


def execute(case):
    work=None;sid=case['source_id'];year=int(case['year']);pol=case['pollutant']
    status=dict(case)
    try:
        days,mass,info,row=emissions(case);annual=math.fsum(mass)
        if annual<CFG['minimum_tonnes']:return dict(status,state='below_threshold',annual_emissions_tonnes=annual)
        coord=coordinates(case,row)
        if coord is None:return dict(status,state='outside_boundary')
        stack,filled,assumption=stack_for(case)
        dists=distances(coord['latitude'],coord['longitude']);gid=ident(GRID[int(np.argmin(dists))]['grid_id'])
        context=met_context(year,gid)
        group=Path(CFG['output_root'])/FOLDERS[CFG['category']]/str(year)/pol;group.mkdir(parents=True,exist_ok=True)
        target=group/f'aermod_{sid}.csv.gz'
        config_values={k:CFG[k] for k in ['version','script_sha256','helper_sha256','executable_sha256','powerplant_mass_unit','stack_file_sha256','elevation_file_sha256','minimum_tonnes','grid_file_sha256','boundary_sha256','aermap_root']}
        fingerprint,legacy_fingerprint=radius_fingerprints(dict(config=config_values,case=case,inventory=info,coordinates=coord,stack=stack,met=context))
        radius=radius_km(annual)
        old_radius=30 if annual>=3600 else 20 if annual>=1000 else 15
        old=core.previous(group,sid)
        if old:
            core.require(old['fingerprint'] in {fingerprint,legacy_fingerprint} and target.is_file() and core.sha(target)==old['output_files'][0]['sha256'],'Existing scratch case differs beyond the radius correction; share the error')
            core.require(all(core.file_info(r['path'])==r for r in old['receptor_inputs']),'Receptor files changed since the completed scratch run')
            if old['fingerprint']==fingerprint or radius==old_radius:
                interrupted=group/(sid+'.running')
                if interrupted.is_dir():core.safe_remove_work(interrupted)
                return dict(status,state='skipped_verified_complete')
        else:core.require(not target.exists(),'Unregistered scratch concentration file; inspect receipt')
        work=group/(sid+'.running')
        if work.exists():
            archive=group/'unfinished_attempts';archive.mkdir(exist_ok=True)
            core.require(not work.is_symlink(),'Unexpected work symlink')
            work.rename(archive/(sid+'_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')))
        work.mkdir()
        points,recfiles=receptor_network(coord['latitude'],coord['longitude'],radius,work)
        release=dict(coord,**stack,source_type='POINT',release_id='S1',weight=1.)
        mass_audit=core.write_hourly(work/'hourly_emissions.txt',days,mass,[release])
        text=control(case,release,context);(work/'aermod.inp').write_text(text,encoding='ascii')
        with (work/'run.log').open('wb') as log:
            result=subprocess.run([CFG['executable'],'aermod.inp','model.out'],cwd=work,stdout=log,stderr=subprocess.STDOUT,timeout=21600,check=False)
        errors=(work/'model.err').read_text(errors='replace') if (work/'model.err').is_file() else ''
        model=core.model_summary(work/'model.out',result.returncode,year,context,errors)
        saved=work/target.name;core.daily_csv(work/'model.out',saved,year,points)
        record=dict(status,state='completed',fingerprint=fingerprint,model=model,met_context=context,inventory=info,
            annual_emissions_tonnes=annual,unique_receptors=len(points),source_release_parameters=[release],
            stack_parameters_filled_from_category_medians=filled,release_type_assumption=assumption,receptor_inputs=recfiles,
            hourly_mass_balance=mass_audit,output_files=[dict(path=str(target),size=saved.stat().st_size,sha256=core.sha(saved))])
        record['receptor_radius']=dict(km=radius,annual_short_tons=annual/.90718474,threshold_basis='US short tons/year')
        if old:
            backup=group/'previous_radius_results'/target.name
            record.update(previous_radius_fingerprint=old['fingerprint'],previous_radius_output=old['output_files'][0],previous_radius_csv=str(backup))
            core.require(not backup.exists(),'Previous-radius backup already exists; inspect before replacement')
        (work/'receipt.json').write_text(json.dumps(record))
        if old:
            backup.parent.mkdir(parents=True,exist_ok=True);target.rename(backup)
        saved.rename(target);commit_result(group,sid,record,work)
        core.safe_remove_work(work)
        return dict(status,state='completed',annual_emissions_tonnes=annual,filled_parameters=filled)
    except Exception as exc:
        if work and work.is_dir():(work/'failure.txt').write_text(traceback.format_exc())
        return dict(status,state='failed',error=str(exc),work_folder=str(work) if work else None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--category',choices=['power_plants'],default='power_plants')
    parser.add_argument('--workers',type=int,default=40)
    parser.add_argument('--data-root',type=Path,default=BASE/'.')
    parser.add_argument('--output-root',type=Path,default=BASE/'outputs/background_model_runs')
    parser.add_argument('--powerplant-mass-unit',choices=['short_ton','metric_tonne'],default='short_ton')
    parser.add_argument('--minimum-tonnes',type=float,default=0.)
    parser.add_argument('--years',type=int,nargs='+',default=list(range(2004,2024)))
    parser.add_argument('--pollutants',nargs='+',choices=POLS,default=POLS)
    parser.add_argument('--stack-file',type=Path)
    parser.add_argument('--elevation-file',type=Path)
    parser.add_argument('--prepare-only',action='store_true',help='Write the input review and case plan without running AERMOD')
    args=parser.parse_args();power=args.category=='power_plants'
    core.require(1<=args.workers<=int(os.environ.get('SLURM_CPUS_PER_TASK','1')),'Workers exceed allocated CPUs')
    core.require(math.isfinite(args.minimum_tonnes) and args.minimum_tonnes>=0,'Invalid minimum mass')
    stack=args.stack_file or (BASE/'inputs/parameters/background/power_plant_stack_parameters.csv' if power else BASE/'inputs/parameters/background/industrial_stack_parameters.csv')
    elev=args.elevation_file or (BASE/'inputs/parameters/background/power_plant_elevation.csv' if power else BASE/'inputs/parameters/background/industrial_elevation.csv')
    output=args.output_root.resolve();data=args.data_root.resolve()
    core.require(not data.is_relative_to(output),'Use scratch outside the group directory for new background results')
    output.mkdir(parents=True,exist_ok=True)
    import fcntl
    lock=(output/(args.category+'.lock')).open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise RuntimeError('This background category is already running')
    grid=BASE/'inputs/supporting/grids/meteorological_grid_10km.csv'
    boundary=data/'inputs/supporting/boundary/study_boundary.shp'
    exe=BASE/'executables/aermod/aermod'
    config=dict(version=VERSION,category=args.category,data_root=str(data),output_root=str(output),
        stack_file=str(stack),elevation_file=str(elev),stack_file_sha256=core.sha(stack),elevation_file_sha256=core.sha(elev),
        met_root=str(BASE/'inputs/meteorology/aermet'),aermap_root=str(BASE/'inputs/supporting/receptors'),
        dem_root=str(BASE/'inputs/supporting/elevation'),grid_file=str(grid),boundary=str(boundary),
        grid_file_sha256=core.sha(grid),boundary_sha256={s:core.sha(boundary.with_suffix(s)) for s in ['.shp','.shx','.dbf','.prj']},
        executable=str(exe),executable_sha256=core.sha(exe),script_sha256=core.sha(Path(__file__)),
        helper_sha256=core.sha(Path(core.__file__)),powerplant_mass_unit=args.powerplant_mass_unit,minimum_tonnes=args.minimum_tonnes)
    initialize(config);tasks=[];issues=[];seen=set();omitted=0;omitted_records=[]
    # Build from the current inventory, not old modeled folders. Newly added
    # non-UOGD facilities are therefore included when their inputs are complete.
    for year in args.years:
        for pol in args.pollutants:
            path=inventory_path(year,pol);table,_=inventory(str(path))
            for facility,unit in table:
                sid=facility+'_'+unit if power else facility
                core.require(re.fullmatch(r'[A-Za-z0-9_-]+',sid),'Unsafe facility/unit identifier: '+sid)
                case=dict(source_category=args.category,source_id=sid,facility_id=facility,unit_id=unit,year=str(year),pollutant=pol)
                try:
                    _,mass,_,row=emissions(case);annual=math.fsum(mass)
                    if annual<=0 or annual<args.minimum_tonnes:
                        omitted+=1;omitted_records.append(dict(case,reason='zero_mass' if annual<=0 else 'below_selected_threshold',annual_emissions_tonnes=annual));continue
                    coord=coordinates(case,row)
                    if coord is None:
                        omitted+=1;omitted_records.append(dict(case,reason='outside_study_boundary',annual_emissions_tonnes=annual));continue
                    stack_for(case)
                    key=(year,pol,sid);core.require(key not in seen,'Repeated case');seen.add(key);tasks.append(case)
                except Exception as exc:issues.append(dict(case,issue=str(exc)))
            print(f'Prepared {year} {pol}: selected={len(tasks):,}; input issues={len(issues):,}',flush=True)
    report=output/'preparation'/args.category;report.mkdir(parents=True,exist_ok=True)
    for name,records,cols in [('input_issues.csv',issues,['source_category','source_id','facility_id','unit_id','year','pollutant','issue']),
                              ('omitted_cases.csv',omitted_records,['source_category','source_id','facility_id','unit_id','year','pollutant','reason','annual_emissions_tonnes']),
                              ('case_plan.csv',tasks,['source_category','source_id','facility_id','unit_id','year','pollutant'])]:
        with (report/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(records)
    (report/'SETUP.json').write_text(json.dumps(dict(config,selected=len(tasks),input_issues=len(issues),zero_below_threshold_or_outside=omitted),indent=2))
    if args.prepare_only:
        print('Preparation only; no AERMOD runs started. Review:',report,flush=True)
        return
    core.require(not issues,'Input preparation needs review before modeling; share '+str(report/'input_issues.csv'))
    core.require(tasks,'No eligible background cases')
    logdir=output/'logs'/args.category;logdir.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f');logpath=logdir/(stamp+'.jsonl');summarypath=logdir/(stamp+'_SUMMARY.json')
    counts={};finished=0;submitted=0;queue=iter(tasks);exhausted=False
    print(f'Selected {len(tasks):,} cases; workers={args.workers}; OUTPUT={output}',flush=True)
    with logpath.open('x') as log,ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,initargs=(config,)) as pool:
        pending={}
        while True:
            while not exhausted and counts.get('failed',0)<20 and len(pending)<args.workers*2:
                try:case=next(queue)
                except StopIteration:exhausted=True;break
                group=output/FOLDERS[args.category]/case['year']/case['pollutant']
                if group.is_dir():recover_radius_update(group,case['source_id'])
                pending[pool.submit(execute,case)]=case;submitted+=1
            if not pending:break
            done,_=wait(pending,return_when=FIRST_COMPLETED)
            for future in done:
                case=pending.pop(future)
                try:record=future.result()
                except Exception as exc:record=dict(case,state='failed',error=str(exc))
                counts[record['state']]=counts.get(record['state'],0)+1;finished+=1
                log.write(json.dumps(record)+'\n');log.flush()
                if record['state']=='failed':print('FAILED',case['source_id'],record['error'],flush=True)
                if finished%100==0:print(f'{finished:,}/{len(tasks):,}: {counts}',flush=True)
                summary=dict(category=args.category,selected=len(tasks),finished=finished,submitted=submitted,counts=counts,
                    all_complete=finished==len(tasks) and not counts.get('failed'),log=str(logpath),output_root=str(output))
                tmp=summarypath.with_suffix('.tmp');tmp.write_text(json.dumps(summary,indent=2));tmp.replace(summarypath)
    print('SUMMARY:',summarypath,flush=True)
    core.require(finished==len(tasks) and not counts.get('failed'),'Some background cases remain incomplete; inspect log')


if __name__=='__main__':main()
