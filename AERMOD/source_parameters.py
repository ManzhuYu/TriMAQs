"""Apply model-ready source parameters to each AERMOD simulation."""
from pathlib import Path
import json,hashlib,re,copy,csv,calendar,math
from datetime import datetime,timedelta

ENGINE=None
ACTIVE=None
PARAMETERS=None
PARAMETER_SHA=None


def initialize_parameters():
    global PARAMETERS,PARAMETER_SHA
    if PARAMETERS is not None:return
    from trimaqs_config import path
    folder=Path(path('source_parameter_folder'))
    site_files={'drilling_completion':'drilling_completion_data.csv','production':'production_data.csv',
                'compressors':'compressor_data_with_elevation.csv','other_tier2':'other_tier2_data_with_elevation.csv'}
    names=list(site_files.values())+['stack_parameters.csv','model_cases.csv']
    digest=hashlib.sha256()
    for name in names:
        with (folder/name).open('rb') as handle:
            for block in iter(lambda:handle.read(1024*1024),b''):digest.update(block)
    PARAMETER_SHA=digest.hexdigest();locations={};sets={};owners={};cases={};parameters={}
    for category,name in site_files.items():
        with (folder/name).open(newline='',encoding='utf-8-sig') as handle:
            for row in csv.DictReader(handle):
                key=(category,row['location_id'])
                if key in locations:raise ValueError('Duplicate location key')
                locations[key]={k:float(row[k]) for k in ['latitude','longitude','elevation_m']}
                locations[key]['location_note']=row.get('location_note','')
                locations[key]['source_id']=row['source_id']
    with (folder/'stack_parameters.csv').open(newline='',encoding='utf-8-sig') as handle:
        for row in csv.DictReader(handle):
            release={k:float(row[k]) for k in ['weight','latitude','longitude','easting_m','northing_m','elevation_m','height_m','temperature_K','velocity_m_s','diameter_m','x_length_m','y_length_m','angle_degrees'] if row.get(k) not in {'',None}}
            release.update(model_source_id=row['release_id'],source_type=row['source_type'],reference_year=row.get('reference_year',''),parameter_method=row.get('parameter_method',''),release_id=row.get('reference_release_id',''),reference_release_pollutant=row.get('reference_pollutant',''),SCCs=row.get('scc_codes','').split(';') if row.get('scc_codes') else [])
            owner=(row['source_category'],row['source_id'])
            if row['release_set_id'] in owners and owners[row['release_set_id']]!=owner:raise ValueError('Conflicting release-set ownership')
            owners[row['release_set_id']]=owner
            if any(r['model_source_id']==release['model_source_id'] for r in sets.get(row['release_set_id'],[])):raise ValueError('Repeated release identifier')
            sets.setdefault(row['release_set_id'],[]).append(release)
    with (folder/'model_cases.csv').open(newline='',encoding='utf-8-sig') as handle:
        for row in csv.DictReader(handle):
            category=row['source_category'];workflow='activity' if category=='drilling_completion' else category
            key='|'.join([workflow,row['year'],row['pollutant'],row['source_id']])
            if key in cases:raise ValueError('Duplicate model case')
            location=locations[(category,row['location_id'])]
            if location['source_id']!=row['source_id'] or owners[row['release_set_id']]!=(category,row['source_id']):raise ValueError('Source/location/release relationship differs')
            entry=dict(location,releases=sets[row['release_set_id']],
                met_grid_id=row['meteorological_grid_id'],receptor_grid_ids=row['receptor_grid_ids'].split(';'),
                annual_emissions_tonnes=float(row['annual_emissions_tonnes']))
            cases[key]=key;parameters[key]=entry
    PARAMETERS={'cases':cases,'parameters':parameters}


def case_parameters():
    initialize_parameters()
    year,pollutant,source=ACTIVE
    workflow=ENGINE['WORKFLOW']
    key='|'.join([workflow,str(year),pollutant,source])
    if key not in PARAMETERS['cases']:raise ValueError('Model parameters missing: '+key)
    return PARAMETERS['parameters'][PARAMETERS['cases'][key]]


def coordinates(source,row):
    item=case_parameters()
    return item['latitude'],item['longitude'],item['elevation_m']


def stack_values(source=None):
    r=case_parameters()['releases'][0]
    return dict(stack_height_m=r['height_m'],stack_diameter_m=r.get('diameter_m',1.),
                exit_temp_k=r.get('temperature_K',300.),exit_velocity_m_s=r.get('velocity_m_s',1.))


def fixed_stack():
    s=stack_values()
    return tuple(s[k] for k in ['stack_height_m','stack_diameter_m','exit_temp_k','exit_velocity_m_s'])


def fingerprint(task):
    initialize_parameters();entry=case_parameters()
    return dict(task=list(task),inventory=ENGINE['file_record'](ENGINE['inventory_path'](task[0],task[1])),
                parameters_sha256=PARAMETER_SHA,source_parameters=entry,
                script_sha256=hashlib.sha256(Path(ENGINE['__file__']).read_bytes()).hexdigest(),
                source_writer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                mass_grams_per_ton=ENGINE['MASS_GRAMS_PER_TON'],
                supporting_inputs=[ENGINE['file_record'](ENGINE[k]) for k in ['AERMOD_EXECUTABLE','GRID_CSV_FILE']])


def write_hourly(folder,source,year,daily,stack):
    releases=case_parameters()['releases']
    if abs(sum(r['weight'] for r in releases)-1)>1e-10:raise ValueError('Release weights do not sum to one')
    filename=folder/'hourly_emissions.txt';expected=0.;written=0.
    with filename.open('w',encoding='ascii') as handle:
        for day,rate in sorted(daily.items()):
            yy,mm,dd=map(int,day.split('-'));expected+=rate*86400
            records=[]
            for r in releases:
                area=r['x_length_m']*r['y_length_m'] if r['source_type']=='AREA' else 1.
                if area<=0:raise ValueError('Nonpositive release area')
                q=f"{rate*r['weight']/area:.12E}";written+=float(q)*area*86400
                extra=f" {r['temperature_K']:.8f} {r['velocity_m_s']:.8f}" if r['source_type']=='POINT' else ''
                records.append((r['model_source_id'],q,extra))
            for hour in range(1,25):
                for sid,q,extra in records:handle.write(f'SO HOUREMIS {yy:04d} {mm:02d} {dd:02d} {hour:02d} {sid} {q}{extra}\n')
    if abs(written-expected)>max(1e-6,abs(expected)*1e-9):raise ValueError('Hourly mass balance failed')
    return str(filename)


def write_input(*args):
    filename=ENGINE['_write_control_file'](*args)
    if not filename:raise ValueError('AERMOD input writer failed')
    path=Path(filename);text=path.read_text()
    block=re.search(r'(?ms)^SO STARTING\s*\n.*?^SO FINISHED',text)
    if block is None:raise ValueError('AERMOD source section not found')
    hourly=args[2] if ENGINE['WORKFLOW'] in {'activity','production'} else args[3]
    extras=[]
    for line in block[0].splitlines():
        fields=line.split()
        if fields and fields[0] in {'GASDEPOS','PARTDIAM','PARTDENS','MASSFRAX'}:extras.append(fields)
    lines=['SO STARTING','   ELEVUNIT METERS']
    for r in case_parameters()['releases']:
        sid=r['model_source_id']
        lines.append(f"   LOCATION {sid} {r['source_type']} {r['easting_m']:.8f} {r['northing_m']:.8f} {r['elevation_m']:.8f}")
        if r['source_type']=='POINT':
            values=[r[k] for k in ['height_m','temperature_K','velocity_m_s','diameter_m']]
        elif r['source_type']=='AREA':
            values=[r[k] for k in ['height_m','x_length_m','y_length_m','angle_degrees']]
        else:raise ValueError('Unsupported release type')
        lines.append('   SRCPARAM '+sid+' 0.0 '+' '.join(f'{v:.8f}' for v in values))
        lines.append(f'   HOUREMIS "{hourly}" {sid}')
        for fields in extras:lines.append('   '+' '.join([fields[0],sid,*fields[2:]]))
    lines+=['   SRCGROUP ALL','SO FINISHED']
    path.write_text(text[:block.start()]+'\n'.join(lines)+text[block.end():],encoding='ascii')
    return filename


def meteorology(folder,gid,year):
    return ENGINE['_prepare_meteorology'](folder,case_parameters()['met_grid_id'],year)


def receptors(folder,lat,lon,radius):
    points={};files=[]
    for gid in case_parameters()['receptor_grid_ids']:
        filename=Path(ENGINE['AERMAP_GRIDWISE_DIR'])/f'Grid_{gid}'/f'receptor_output_{gid}_grid.DAT'
        files.append(ENGINE['file_record'](filename))
        for line in filename.read_text().splitlines():
            fields=line.split()
            if fields and fields[0]=='RE':fields=fields[1:]
            if not fields or fields[0] not in {'EVALCART','DISCCART'}:continue
            xy=tuple(round(float(v),2) for v in fields[1:3])
            if xy in points and points[xy].split()[1:6]!=fields[1:6]:raise ValueError('Conflicting repeated receptor heights')
            points.setdefault(xy,' '.join(fields))
    if not points:raise ValueError('No current receptor points')
    filename=folder/'receptors.DAT'
    filename.write_text('RE ELEVUNIT METERS\n'+'\n'.join('RE '+line for line in points.values())+'\n',encoding='ascii')
    return str(filename),set(points),files


def run_case(task):
    global ACTIVE
    ACTIVE=task
    # Published inventories contain only sources selected for the study.
    result=ENGINE['_run_model_case'](task)
    if result.get('state')=='completed':
        folder=Path(result['folder']);filename=folder/'AERMOD_DONE.json'
        record=json.loads(filename.read_text());record.pop('stack_parameters',None)
        record['source_release_parameters']=copy.deepcopy(case_parameters()['releases'])
        ENGINE['dump'](filename,record)
    return result


def inventory(year,pollutant):
    p='NOX' if pollutant in {'NOX','NO2'} else pollutant
    token='pm25' if p=='PM2.5' else p.lower()
    return ENGINE['INVENTORY_TEMPLATE'].format(year=year,pollutant=p,pollutant_slug=token)


def id_column():return 'source_id'


def emissions_for(row,year):
    days=366 if calendar.isleap(year) else 365
    daily={};annual=0.
    if ENGINE['WORKFLOW']=='activity':
        for i in range(days):
            label=(datetime(year,1,1)+timedelta(days=i)).strftime('%Y-%m-%d')
            value=ENGINE['numeric'](row.get(label))
            if value<0:raise ValueError('Negative daily mass')
            daily[label]=value*ENGINE['MASS_GRAMS_PER_TON']/86400;annual+=value
    else:
        for month in range(1,13):
            label=calendar.month_abbr[month].lower()+'_tonnes';value=ENGINE['numeric'](row.get(label))
            if value<0:raise ValueError('Negative monthly mass')
            length=calendar.monthrange(year,month)[1];annual+=value
            for day in range(1,length+1):daily[f'{year}-{month:02d}-{day:02d}']=value*ENGINE['MASS_GRAMS_PER_TON']/(86400*length)
    if annual>0:
        expected=case_parameters().get('annual_emissions_tonnes',annual)
        if not math.isclose(annual,expected,rel_tol=1e-9,abs_tol=0.):raise ValueError('Inventory annual mass differs from the model case table')
    radius=(20 if annual>=3600 else 15 if annual>=1000 else 12) if ENGINE['WORKFLOW']=='production' else (50 if annual>=3600 else 35 if annual>=1000 else 25)
    return daily,annual,float(radius),0


def configure(namespace):
    """Connect the common source reader to the selected category driver."""
    global ENGINE
    ENGINE=namespace
