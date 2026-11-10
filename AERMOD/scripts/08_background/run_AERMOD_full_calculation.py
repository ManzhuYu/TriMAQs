"""Run the four report pollutants without FASTALL, writing only to scratch.

Uses the current shared inventory/case/release tables and saved physical
control files. Successful cases keep one CSV.gz; one SQLite database per
category/year/pollutant holds completion records, inputs and diagnostics.
Run as a saved Python file with the supplied Slurm jobs. Requires numpy.
"""
from pathlib import Path
import os
for _name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[_name]='1'
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from datetime import datetime,timedelta
from functools import lru_cache
from contextlib import closing,contextmanager
import argparse,calendar,csv,gzip,hashlib,json,math,re,shlex,shutil,sqlite3,subprocess,traceback,zlib
import numpy as np

VERSION='2026-11-07-full-calculation-v1'
BASE=Path(os.environ.get('TRIMAQS_DATA_ROOT',str(Path(__file__).resolve().parents[3]/'TriMAQs_data')))
FOLDERS={'drilling_completion':'Well_drilling_new_4','production':'Wellpad_production3',
         'compressors':'Compressors_4','other_tier2':'ONG_related_4'}
PREFIX={'drilling_completion':'drilling_completion','production':'production',
        'compressors':'compressor','other_tier2':'other_tier2'}
REPORT_POLLUTANTS=['NO2','PM2.5','SO2','VOC']
CONFIG={}
RELEASES={}


def require(ok,message):
    if not ok:raise ValueError(message)


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def signature(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def file_info(path):
    path=Path(path).resolve();stat=path.stat()
    return dict(path=str(path),size=stat.st_size,mtime_ns=stat.st_mtime_ns)


def fields(line):
    p=shlex.split(line)
    return p[1:] if p and p[0] in {'CO','SO','RE','ME','OU'} else p


def initialize(config):
    global CONFIG,RELEASES
    CONFIG=config;RELEASES={}
    path=Path(config['data_root'])/'inputs/parameters/stack_parameters.csv'
    with path.open(newline='',encoding='utf-8-sig') as handle:
        for row in csv.DictReader(handle):
            if row['source_category']==config['category']:
                RELEASES.setdefault(row['release_set_id'],[]).append(row)


def inventory_path(case):
    pol='NOX' if case['pollutant']=='NO2' else case['pollutant']
    token='pm25' if pol=='PM2.5' else pol.lower()
    return Path(CONFIG['data_root'])/'inputs/emissions'/case['source_category']/pol/f"{PREFIX[case['source_category']]}_emissions_{token}_{case['year']}.csv"


@lru_cache(maxsize=4)
def inventory_table(name):
    path=Path(name);before=file_info(path);digest=sha(path);rows={}
    with path.open(newline='',encoding='utf-8-sig') as handle:
        for row in csv.DictReader(handle):
            sid=row['source_id'];require(sid not in rows,'Repeated inventory source: '+sid);rows[sid]=row
    require(file_info(path)==before,'Inventory changed while reading: '+str(path))
    return rows,dict(before,sha256=digest)


def daily_mass(case):
    table,info=inventory_table(str(inventory_path(case)))
    require(case['source_id'] in table,'Source absent from current shared inventory')
    row=table[case['source_id']];year=int(case['year'])
    n=366 if calendar.isleap(year) else 365
    days=[datetime(year,1,1)+timedelta(days=i) for i in range(n)]
    if case['source_category']=='drilling_completion':
        mass=[float(row[d.strftime('%Y-%m-%d')]) for d in days]
    else:
        mass=[float(row[calendar.month_abbr[d.month].lower()+'_tonnes'])/calendar.monthrange(year,d.month)[1] for d in days]
    require(all(math.isfinite(v) and v>=0 for v in mass),'Invalid shared inventory mass')
    require(math.isclose(math.fsum(mass),float(case['annual_emissions_tonnes']),rel_tol=1e-9,abs_tol=1e-12),
            'Inventory annual mass differs from current model_cases.csv')
    require(math.fsum(mass)>0,'Nonpositive selected case')
    return days,mass,info


def resolve_reference(folder,value):
    path=Path(value)
    if not path.is_absolute():path=folder/path
    if path.is_file():return path.resolve()
    local=folder/path.name
    require(local.is_file(),'Saved physical input missing: '+str(path))
    return local.resolve()


@lru_cache(maxsize=4)
def baseline_records(group):
    db=Path(group)/'audit.sqlite3'
    if not db.is_file():return {}
    with closing(sqlite3.connect(db.resolve().as_uri()+'?mode=ro',uri=True,timeout=120)) as conn:
        return {sid:json.loads(record) for sid,record in conn.execute('SELECT source,record FROM cases')}


def baseline_record(folder):
    direct=folder/'AERMOD_DONE.json'
    record=json.loads(direct.read_text()) if direct.is_file() else baseline_records(str(folder.parent)).get(folder.name)
    require(record and record.get('state')=='completed','Original completed case/audit missing: '+str(folder))
    return record


def physical_inputs(case,releases):
    folder=Path(CONFIG['case_root'])/FOLDERS[case['source_category']]/case['year']/case['pollutant']/case['source_id']
    files=list(folder.glob('*.inp'));require(len(files)==1,'Expected one saved model input: '+str(folder))
    text=files[0].read_text(encoding='ascii');options=[];pollutants=[];locations={};parameters={};refs={};includes=[];periods=[]
    block=re.search(r'(?ms)^RE STARTING\s*\n.*?^RE FINISHED',text)
    require(block is not None,'Saved receptor block absent')
    for line in text.splitlines():
        p=fields(line)
        if not p:continue
        if p[0]=='MODELOPT':options.append(p[1:])
        if p[0]=='POLLUTID':pollutants.append(p[1])
        if p[0]=='LOCATION':locations[p[1]]=p[2:]
        if p[0]=='SRCPARAM':parameters[p[1]]=p[2:]
        if p[0]=='STARTEND':periods.append(p[1:])
        if p[0] in {'SURFFILE','PROFFILE'}:
            require(p[0] not in refs,'Repeated meteorological input');refs[p[0]]=resolve_reference(folder,p[1])
    require(len(options)==1 and 'FASTALL' in options[0] and pollutants==[case['pollutant']],
            'Saved pollutant/options differ from the expected FASTALL baseline')
    year=int(case['year'])
    require(periods==[[str(year),'01','01',str(year),'12','31']], 'Expected a full annual saved modeling period')
    require(set(refs)=={'SURFFILE','PROFFILE'},'Meteorological input pair absent')
    require(set(locations)=={r['release_id'] for r in releases}==set(parameters),'Shared release IDs differ from saved model')
    for r in releases:
        sid=r['release_id'];typ=r['source_type'];loc=locations[sid];par=parameters[sid]
        require(typ in {'POINT','AREA'} and loc[0]==typ,'Unsupported or changed source type')
        keys=['height_m','temperature_K','velocity_m_s','diameter_m'] if typ=='POINT' else ['height_m','x_length_m','y_length_m','angle_degrees']
        require(len(loc)>=4 and len(par)>=5,'Incomplete saved release')
        for actual,key in zip(loc[1:4],['easting_m','northing_m','elevation_m']):
            require(math.isclose(float(actual),float(r[key]),rel_tol=1e-8,abs_tol=1e-5),'Changed shared source geometry: '+key)
        for actual,key in zip(par[1:5],keys):
            require(math.isclose(float(actual),float(r[key]),rel_tol=1e-8,abs_tol=1e-5),'Changed shared stack parameter: '+key)
    for line in block[0].splitlines():
        p=fields(line)
        if p and p[0]=='INCLUDED':includes.append(resolve_reference(folder,p[1]))
    require(includes,'Expected saved receptor include file')
    points=set()
    for path in includes:
        for line in path.read_text().splitlines():
            p=fields(line)
            require(not p or p[0]!='INCLUDED','Nested receptor include not supported')
            if p and p[0] in {'EVALCART','DISCCART'}:points.add((round(float(p[1]),2),round(float(p[2]),2)))
    require(points,'No saved Cartesian receptors')
    old=baseline_record(folder);context=old['met_context']
    require(str(context['grid_id'])==str(case['meteorological_grid_id']) and int(context['year'])==year,'Saved meteorological grid/year differs')
    require(len(points)==old['unique_receptors'],'Saved receptor count differs from completed case')
    for key,name in [('SURFFILE','output_stage2.SFC'),('PROFFILE','output_stage2.PFL')]:
        info=file_info(refs[key]);expected=context['files'][name]
        require(info['size']==expected['size'] and info['mtime_ns']==expected['mtime_ns'],'Meteorology changed since completed baseline: '+name)
    return text,files[0],refs,includes,points,context


def write_hourly(path,days,mass,releases):
    weights=[float(r['weight']) for r in releases]
    require(all(math.isfinite(w) and w>=0 for w in weights) and abs(math.fsum(weights)-1)<1e-10,'Invalid release allocation weights')
    written=0.
    with path.open('x',encoding='ascii') as output:
        for day,total in zip(days,mass):
            records=[]
            for r,weight in zip(releases,weights):
                area=float(r['x_length_m'])*float(r['y_length_m']) if r['source_type']=='AREA' else 1.
                require(math.isfinite(area) and area>0,'Invalid area release dimensions')
                q=f'{total*1e6/86400*weight/area:.12E}'
                written+=float(q)*area*86400/1e6
                extra=f" {float(r['temperature_K']):.8f} {float(r['velocity_m_s']):.8f}" if r['source_type']=='POINT' else ''
                records.append((r['release_id'],q,extra))
            for hour in range(1,25):
                for sid,q,extra in records:output.write(f'SO HOUREMIS {day:%Y %m %d} {hour:02d} {sid} {q}{extra}\n')
    require(math.isclose(written,math.fsum(mass),rel_tol=1e-9,abs_tol=1e-12),'Hourly emission mass balance failed')
    return dict(inventory_tonnes=math.fsum(mass),hourly_tonnes=written,
                daily_mass_sha256=signature(mass),allocation='Shared release weights; uniform within each inventory day')


def control_text(original,work,refs,includes):
    output=[];include_index=0;count=0
    for line in original.splitlines():
        p=fields(line)
        if not p:output.append(line);continue
        key=p[0]
        if key in {'ERRORFIL','SUMMFILE'}:continue
        if key=='MODELOPT':
            require('FASTALL' in p[1:],'Saved FASTALL option absent')
            line='   MODELOPT '+' '.join(x for x in p[1:] if x!='FASTALL');count+=1
        elif key in refs:line='   '+key+' "'+str(refs[key])+'"'
        elif key=='HOUREMIS':line='   HOUREMIS "hourly_emissions.txt" '+' '.join(p[2:])
        elif key=='INCLUDED':
            require(include_index<len(includes),'Unexpected include reference')
            line='   INCLUDED "'+str(includes[include_index])+'"';include_index+=1
        output.append(line)
    require(count==1 and include_index==len(includes),'Ambiguous control/input blocks')
    text='\n'.join(output)+'\n'
    text,n=re.subn(r'(?ms)^OU STARTING\s*\n.*?^OU FINISHED',
        'OU STARTING\n   SUMMFILE "model.sum"\n   DAYTABLE 24\nOU FINISHED',text)
    require(n==1,'Output block absent/ambiguous')
    # Full annual PM2.5 runs retain their original POLLUTID processing.
    # H1H was needed only for the earlier short diagnostic.
    return text.replace('CO FINISHED','   ERRORFIL "model.err"\nCO FINISHED')


def model_summary(path,returncode,year,context,errors):
    fatals=[];missing=[];hours=[];versions=[];summary=[];units=[]
    with path.open(errors='replace') as handle:
        for line in handle:
            for pattern,target in [(r'A Total of\s+(\d+)\s+Fatal Error',fatals),(r'A Total of\s+(\d+)\s+Missing Hours',missing),
                                   (r'A Total of\s+(\d+)\s+Hours Were Processed',hours)]:
                m=re.search(pattern,line,re.I)
                if m:target.append(int(m[1]))
            m=re.search(r'AERMOD\s*-\s*VERSION\s+(\d+)',line,re.I)
            if m:versions.append(m[1])
            if 'A Total of' in line:summary.append(line.strip())
            if 'Output Units' in line:units.append(line.strip())
    messages=[line.strip() for line in errors.splitlines() if re.search(r'^\s*[A-Z]{2}\s+E\d{3}\b',line)]
    require(returncode==0 and fatals and max(fatals)==0 and not messages,'AERMOD failed: '+('; '.join(messages) or str(summary[-15:])))
    require(versions and versions[-1]=='24142','Unexpected AERMOD version')
    require(hours and hours[-1]==24*(366 if calendar.isleap(year) else 365),'Incomplete annual model hours')
    expected=len(context['exceptions'])
    require(missing and missing[-1]==expected,'Missing meteorological hours differ from the previously reviewed exceptions')
    return dict(version=versions[-1],fatal_errors=max(fatals),processed_hours=hours[-1],missing_met_hours_reported=missing[-1],
                units=units,summary_lines=summary[-20:],documented_exceptions=context['exceptions'])


def daily_csv(raw,destination,year,points):
    ordered=sorted(points);index={xy:i for i,xy in enumerate(ordered)}
    n=366 if calendar.isleap(year) else 365
    values=np.full((len(ordered),n),np.nan,dtype=np.float64)
    current=None;group=None;recording=False
    with raw.open(errors='replace') as handle:
        for line in handle:
            if 'CONCURRENT' in line:
                current=None;group=None;recording=False
                m=re.search(r'CONCURRENT\s+24\s*(?:-?HR)?\s+(?:AVERAGE\s+)?CONC(?:ENTRATION)?\s+VALUES.*HOUR\s+(\d+)\s+FOR DAY\s+(\d+)\s+OF\s+(\d{4})',line)
                if m and int(m[1])==24 and int(m[3])==year:
                    current=int(m[2])-1;require(0<=current<n,'Invalid output Julian day')
            if 'FOR SOURCE GROUP:' in line:group=line.split('FOR SOURCE GROUP:',1)[1].strip().split()[0]
            if '*** DISCRETE CARTESIAN RECEPTOR POINTS ***' in line:
                recording=current is not None and group=='ALL';continue
            if '***' in line or '\f' in line:recording=False
            if not recording:continue
            p=line.split()
            if len(p) not in (3,6):continue
            try:nums=[float(v.replace('D','E')) for v in p]
            except ValueError:continue
            for i in range(0,len(nums),3):
                x,y,c=nums[i:i+3];xy=(round(x,2),round(y,2))
                require(xy in index and math.isfinite(c) and c>=0,'Unexpected receptor/concentration')
                row=index[xy];old=values[row,current]
                require(np.isnan(old) or old==c,'Conflicting repeated receptor/day')
                values[row,current]=c
    require(np.isfinite(values).all(),'Incomplete receptor/day output; no values were filled')
    with gzip.open(destination,'xt',newline='',compresslevel=3) as handle:
        writer=csv.writer(handle);writer.writerow(['utm_easting','utm_northing']+[f'day{d}_conc' for d in range(1,n+1)])
        for xy,row in zip(ordered,values):writer.writerow([*xy,*row])
    validate_csv(destination,year,points)


def validate_csv(path,year,points):
    n=366 if calendar.isleap(year) else 365;seen=set()
    with gzip.open(path,'rt',newline='') as handle:
        reader=csv.reader(handle)
        require(next(reader)==['utm_easting','utm_northing']+[f'day{d}_conc' for d in range(1,n+1)],'Saved CSV calendar/schema differs')
        for row in reader:
            require(len(row)==n+2,'Saved CSV row length differs')
            x,y=map(float,row[:2]);xy=(round(x,2),round(y,2));v=np.array(row[2:],dtype=float)
            require(xy in points and xy not in seen and np.isfinite(v).all() and (v>=0).all(),'Saved CSV receptor/value invalid')
            seen.add(xy)
    require(seen==points,'Saved CSV receptor coverage differs')


@contextmanager
def database(group):
    conn=sqlite3.connect(group/'audit.sqlite3',timeout=120)
    conn.execute('CREATE TABLE IF NOT EXISTS cases (source TEXT PRIMARY KEY, record TEXT NOT NULL)')
    conn.execute('CREATE TABLE IF NOT EXISTS diagnostics (source TEXT,name TEXT,content BLOB NOT NULL,PRIMARY KEY(source,name))')
    try:
        with conn:yield conn
    finally:conn.close()


def previous(group,sid):
    with database(group) as conn:
        row=conn.execute('SELECT record FROM cases WHERE source=?',(sid,)).fetchone()
    return json.loads(row[0]) if row else None


def safe_remove_work(work):
    base=Path(CONFIG['output_root']).resolve();resolved=work.resolve()
    require(not work.is_symlink() and resolved.is_relative_to(base) and work.name.endswith('.running'), 'Unexpected temporary cleanup target')
    shutil.rmtree(work)


def execute(case):
    work=None;sid=case['source_id'];year=int(case['year']);pol=case['pollutant']
    status=dict(source_category=case['source_category'],source_id=sid,year=year,pollutant=pol)
    try:
        require(re.fullmatch(r'[A-Za-z0-9_-]+',sid),'Unsafe source identifier')
        releases=RELEASES[case['release_set_id']]
        require(all(r['source_id']==sid for r in releases),'Release ownership differs')
        days,mass,inv=daily_mass(case)
        original,inp,refs,includes,points,context=physical_inputs(case,releases)
        dependencies=[file_info(p) for p in [*refs.values(),*includes]]
        fingerprint=signature(dict(version=VERSION,script_sha256=CONFIG['script_sha256'],executable_sha256=CONFIG['executable_sha256'],
            case=case,releases=releases,inventory=inv,saved_input_sha256=sha(inp),physical_files=dependencies))
        group=Path(CONFIG['output_root'])/FOLDERS[case['source_category']]/str(year)/pol
        group.mkdir(parents=True,exist_ok=True)
        result_file=group/f'aermod_{sid}.csv.gz'
        old=previous(group,sid)
        if old and old['fingerprint']==fingerprint and old['state']=='completed':
            require(result_file.is_file() and result_file.stat().st_size==old['output_files'][0]['size'] and sha(result_file)==old['output_files'][0]['sha256'],
                    'Completed scratch result is missing/changed')
            cleanup=group/(sid+'.running')
            if cleanup.is_dir():safe_remove_work(cleanup)
            return dict(status,state='skipped_verified_complete',folder=str(group))
        require(old is None and not result_file.exists(),'Existing scratch result belongs to different inputs; select a new --output-root')
        work=group/(sid+'.running')
        if work.exists():
            archive=group/'unfinished_attempts';archive.mkdir(exist_ok=True)
            require(not work.is_symlink(),'Unexpected work symlink')
            work.rename(archive/(sid+'_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')))
        work.mkdir()
        mass_audit=write_hourly(work/'hourly_emissions.txt',days,mass,releases)
        new=control_text(original,work,refs,includes)
        (work/'aermod.inp').write_text(new,encoding='ascii')
        with (work/'run.log').open('wb') as log:
            result=subprocess.run([CONFIG['executable'],'aermod.inp','model.out'],cwd=work,
                stdout=log,stderr=subprocess.STDOUT,timeout=CONFIG['timeout'],check=False)
        errors=(work/'model.err').read_text(errors='replace') if (work/'model.err').is_file() else ''
        model=model_summary(work/'model.out',result.returncode,year,context,errors)
        target=work/result_file.name
        daily_csv(work/'model.out',target,year,points)
        # All diagnostics are retained in a single group database, not loose files.
        record=dict(status,state='completed',fingerprint=fingerprint,model=model,met_context=context,
            annual_emissions_tonnes=math.fsum(mass),unique_receptors=len(points),hourly_mass_balance=mass_audit,
            source_release_parameters=releases,inventory=inv,baseline_input=str(inp),model_options='FASTALL removed; all other saved options retained',
            output_files=[dict(path=str(result_file),size=target.stat().st_size,sha256=sha(target))])
        diagnostics={'aermod.inp':new.encode('ascii')}
        for name in ['model.err','model.sum','run.log']:
            if (work/name).is_file():diagnostics[name]=(work/name).read_bytes()
        # Same-filesystem rename. If interrupted between rename and DB commit,
        # the next invocation recovers the small receipt inside the work folder.
        receipt=work/'receipt.json';receipt.write_text(json.dumps(record))
        target.rename(result_file)
        with database(group) as conn:
            conn.execute('INSERT INTO cases VALUES (?,?)',(sid,json.dumps(record)))
            conn.executemany('INSERT INTO diagnostics VALUES (?,?,?)',[(sid,name,zlib.compress(data)) for name,data in diagnostics.items()])
        safe_remove_work(work)
        return dict(status,state='completed',folder=str(group),missing_met_hours=model['missing_met_hours_reported'])
    except Exception as exc:
        if work and work.is_dir():(work/'failure.txt').write_text(traceback.format_exc())
        return dict(status,state='failed',error=str(exc),work_folder=str(work) if work else None)


def recover(group,sid):
    """Complete an interrupted result rename before scheduling the case."""
    work=group/(sid+'.running');receipt=work/'receipt.json';target=group/f'aermod_{sid}.csv.gz'
    if not receipt.is_file() or not target.is_file():return
    record=json.loads(receipt.read_text())
    require(record['source_id']==sid and record['output_files'][0]['path']==str(target),'Interrupted receipt points elsewhere')
    require(target.stat().st_size==record['output_files'][0]['size'] and sha(target)==record['output_files'][0]['sha256'],'Interrupted saved CSV changed')
    with database(group) as conn:
        old=conn.execute('SELECT record FROM cases WHERE source=?',(sid,)).fetchone()
        if old:
            require(json.loads(old[0])['fingerprint']==record['fingerprint'],'Interrupted receipt conflicts with completion DB')
        else:
            conn.execute('INSERT INTO cases VALUES (?,?)',(sid,json.dumps(record)))
            for name in ['aermod.inp','model.err','model.sum','run.log']:
                if (work/name).is_file():conn.execute('INSERT INTO diagnostics VALUES (?,?,?)',(sid,name,zlib.compress((work/name).read_bytes())))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--category',choices=list(FOLDERS),required=True)
    parser.add_argument('--workers',type=int,default=40)
    parser.add_argument('--data-root',type=Path,default=BASE/'.')
    parser.add_argument('--case-root',type=Path,default=BASE/'working/baseline_cases')
    parser.add_argument('--output-root',type=Path,default=BASE/'outputs/model_runs')
    parser.add_argument('--executable',type=Path,default=BASE/'executables/aermod/aermod')
    parser.add_argument('--years',type=int,nargs='+',default=list(range(2004,2024)))
    parser.add_argument('--pollutants',nargs='+',choices=REPORT_POLLUTANTS,default=REPORT_POLLUTANTS)
    parser.add_argument('--timeout',type=int,default=21600)
    parser.add_argument('--stop-after-failures',type=int,default=20)
    args=parser.parse_args()
    require(1<=args.workers<=int(os.environ.get('SLURM_CPUS_PER_TASK','1')),'Workers exceed allocated CPUs; use the batch job')
    require(args.stop_after_failures>0 and args.timeout>0,'Invalid timeout/failure limit')
    output=args.output_root.resolve();source=args.case_root.resolve();data=args.data_root.resolve()
    require(output!=source and not output.is_relative_to(source) and not source.is_relative_to(output),'Scratch output overlaps existing model root')
    require(not output.is_relative_to(data) and not data.is_relative_to(output),'Scratch output overlaps shared dataset')
    require(args.executable.is_file(),'AERMOD executable missing')
    output.mkdir(parents=True,exist_ok=True)
    import fcntl
    lock=(output/(args.category+'.lock')).open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise RuntimeError('This category is already running in the selected scratch output root')
    config=dict(category=args.category,data_root=str(data),case_root=str(source),output_root=str(output),
        executable=str(args.executable.resolve()),executable_sha256=sha(args.executable),
        script_sha256=sha(Path(__file__)),timeout=args.timeout)
    initialize(config)
    with (data/'inputs/parameters/model_cases.csv').open(newline='',encoding='utf-8-sig') as handle:
        cases=[row for row in csv.DictReader(handle) if row['source_category']==args.category and int(row['year']) in args.years and row['pollutant'] in args.pollutants]
    require(cases,'No current shared model cases selected')
    keys=[(r['year'],r['pollutant'],r['source_id']) for r in cases]
    require(len(keys)==len(set(keys)),'Repeated selected model case')
    require(all(r['release_set_id'] in RELEASES for r in cases),'Shared release set absent')
    cases.sort(key=lambda r:(int(r['year']),r['pollutant'],r['source_id']))
    logdir=output/'logs'/args.category;logdir.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    logpath=logdir/(stamp+'.jsonl');summarypath=logdir/(stamp+'_SUMMARY.json')
    print(f'Full calculation: {len(cases):,} selected cases; {args.workers} workers; category={args.category}',flush=True)
    print('OUTPUT:',output,flush=True);print('LOG:',logpath,flush=True)
    counts={'completed':0,'skipped_verified_complete':0,'failed':0};submitted=0;finished=0;queue=iter(cases);exhausted=False
    with logpath.open('x') as log,ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,initargs=(config,)) as pool:
        pending={}
        while True:
            while not exhausted and counts['failed']<args.stop_after_failures and len(pending)<args.workers*2:
                try:case=next(queue)
                except StopIteration:exhausted=True;break
                group=output/FOLDERS[args.category]/case['year']/case['pollutant']
                if group.is_dir():recover(group,case['source_id'])
                pending[pool.submit(execute,case)]=case;submitted+=1
            if not pending:break
            complete,_=wait(pending,return_when=FIRST_COMPLETED)
            for future in complete:
                case=pending.pop(future)
                try:record=future.result()
                except Exception as exc:record=dict(case,state='failed',error=str(exc))
                counts[record['state']]+=1;finished+=1
                log.write(json.dumps(record)+'\n');log.flush()
                if record['state']=='failed':print('FAILED:',case['year'],case['pollutant'],case['source_id'],record['error'],flush=True)
                if finished%100==0:print(f'{finished:,}/{len(cases):,}: {counts}',flush=True)
                summary=dict(version=VERSION,category=args.category,selected=len(cases),submitted=submitted,finished=finished,counts=counts,
                    all_complete=finished==len(cases) and counts['failed']==0,output_root=str(output),log=str(logpath))
                temp=summarypath.with_suffix('.tmp');temp.write_text(json.dumps(summary,indent=2));temp.replace(summarypath)
    print('SUMMARY:',summarypath,flush=True)
    require(finished==len(cases) and counts['failed']==0,'Run incomplete; successful cases are retained. Inspect log and rerun the same job after resolving failures.')


if __name__=='__main__':main()
