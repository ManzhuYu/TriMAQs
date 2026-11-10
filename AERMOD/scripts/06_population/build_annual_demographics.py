"""Build demographic populations consistent with the installed annual population.
Uses the existing census ZIP and population overlap checkpoint. No shared files change.
"""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import os
import re
import sqlite3
import zipfile
import numpy as np
import pandas as pd

# USER SETTINGS --------------------------------------------------------------
DATA_ROOT = Path(os.environ.get('TRIMAQS_DATA_ROOT', str(Path(__file__).resolve().parents[3]/'TriMAQs_data')))
PROJECT_ROOT = DATA_ROOT

CENSUS_ZIP = DATA_ROOT / 'inputs/source_data/census/Census_population.zip'
POPULATION_BUILD = DATA_ROOT / 'inputs/population_demographics/general'
GENERAL = DATA_ROOT / 'inputs/population_demographics/general'
OUTPUT_DIR = DATA_ROOT / 'inputs/population_demographics/demographics'
WORKERS = 20
# Confirmed definition: income-to-poverty ratio below 2.00.
LOW_INCOME_RATIO = 2.0
# ---------------------------------------------------------------------------
VERSION = 'demographic-population-v1'
OUTPUT_POLICY = 'income-availability-v2'
# Source allocations have not changed; retain their existing checkpoint identity.
PREVIOUS_README_HASHES = {
 'README.md': '78d55558494c6b407313ff9ac03624baea51ee055118de447b32986f610b5fea',
 'Low_Income/README.md': '125c9b822bea2a31c103939654ab7357ac8e722f9773f91ef00ac5c42c302c95'
}
AGE = ['male_under5','male_5_64','male_65plus','female_under5','female_5_64','female_65plus']
RACE = ['white_non_hispanic','black_non_hispanic','aian_asian_pacific_non_hispanic','other_non_hispanic','hispanic']
INCOME = ['low_income','other_poverty_determined','poverty_not_determined']
GROUPS = AGE + RACE + INCOME
READMES = {
 'README.md': '''# Annual demographic population

These files contain modeled demographic populations for each 1 km TriMAQs study-grid cell for 2004–2023. The grid IDs and coordinates match the annual total-population and concentration datasets. Counts are in people and are nonnegative whole numbers.

The source data are NHGIS census tract tables for sex and age, race and Hispanic origin, and income-to-poverty ratio. Census 2000 and ACS five-year source populations are allocated using the same source-specific tract–grid overlap weights used for total population. ACS periods are assigned to midpoint years. The incomplete 2005–2009 ACS source is omitted.

For 2004–2007, age/sex and race counts interpolate between the 2000 census and the 2006–2010 ACS estimate assigned to 2008. The supplied archive does not provide a compatible 2000 poverty-ratio anchor for this calculation, so the 2008 grid-level poverty-category shares are held constant for 2004–2007. These early low-income estimates reflect changes in total population, not measured changes in income distribution.

For 2008–2021, categories are allocated from the corresponding ACS release. For 2022–2023, category trends are fitted to the 2017–2021 counts with the intercept fixed at 2021 and slopes restricted to keep fitted values nonnegative over 2017–2023. The resulting category shares are reconciled to the installed population in each grid cell. Total-population assumptions about source periods and recent missing-tract boundaries also apply.

Age and sex are rounded as one joint partition. Race/ethnicity and poverty categories are each rounded as separate partitions. Within every cell, male plus female and the five race/ethnicity categories each equal total population. Under-5 plus age-65-and-over populations do not exceed total population. Poverty categories also sum to total population; people without a determined poverty status remain a separate category. If an early-year cell has residents but its 2008 reference has no population, its income distribution is unavailable. Those residents are retained in poverty_not_determined_population and listed in income_distribution_unavailable.csv; they are not assigned an income category. Rounding preserves each cell's installed population. Subgroup regional totals are modeled estimates, without independent annual demographic controls.

Every CSV contains grid_ID, latitude and longitude, followed by population counts. Group definitions and column names are described in each subfolder. Age, sex, race/ethnicity and income are separate classifications; their counts must not be added across classifications.
''',
 'Age_less_5/README.md': '''# Population under age 5

`Age_less_then_5_YEAR.csv` contains grid_ID, latitude, longitude and under5_population, in people. Counts include both sexes and ages below 5 years. They are derived from the joint age/sex allocation described in the parent README and do not exceed the installed total population. The data cover 2004–2023.
''',
 'Age_more_64/README.md': '''# Population age 65 and over

`Age_over_64_YEAR.csv` contains grid_ID, latitude, longitude and age65_population, in people. Counts include both sexes and ages 65 and over. They are derived from the joint age/sex allocation described in the parent README. Under-5 and age-65-and-over counts together do not exceed the installed total population. The data cover 2004–2023.
''',
 'Gender/README.md': '''# Population by sex

`Sex_YEAR.csv` contains grid_ID, latitude, longitude, male_population and female_population, in people. Categories follow the male/female sex classifications in the source census tables; they do not describe gender identity. The two counts sum to the installed total population in every cell. They are rounded jointly with the age groups, as described in the parent README. The data cover 2004–2023.
''',
 'Race-ethnicity/README.md': '''# Population by race and Hispanic origin

`Race_YEAR.csv` contains grid_ID, latitude and longitude, followed by five mutually exclusive population counts in people:

| Column | Population counted |
|---|---|
| white_non_hispanic | Non-Hispanic White alone |
| black_non_hispanic | Non-Hispanic Black or African American alone |
| aian_asian_pacific_non_hispanic | Non-Hispanic American Indian/Alaska Native, Asian, or Native Hawaiian/other Pacific Islander alone, combined |
| other_non_hispanic | Non-Hispanic some other race alone or two or more races |
| hispanic | Hispanic or Latino, of any race |

The combined category follows the level of detail available consistently in the supplied ACS time-series table. It cannot be interpreted as Asian population alone. Hispanic population is separate from the four non-Hispanic groups, so the five columns sum to total population without double counting. The allocation and temporal methods are described in the parent README. The data cover 2004–2023.
''',
 'Low_Income/README.md': '''# Population below twice the poverty threshold

`Population_low_income_YEAR.csv` contains grid_ID, latitude and longitude, followed by:

| Column | Population counted, in people |
|---|---|
| low_income_population | Income-to-poverty ratio below 2.00 |
| other_poverty_determined_population | Ratio 2.00 or greater |
| poverty_not_determined_population | People outside the source poverty-status universe, plus residents with an unavailable early-year income distribution |
| poverty_determined_population | Sum of the first two counts |

The first three counts sum to installed total population. People whose poverty status is not determined are not classified as low income or as income at or above twice the threshold. If an early-year cell has residents but the 2008 reference cell has no population, no income distribution can be inferred from that reference. Those residents remain unclassified and are listed in income_distribution_unavailable.csv in the parent folder. A zero low-income count in such a cell does not establish an absence of low-income residents. For a low-income share among people with determined poverty status, divide low_income_population by poverty_determined_population and leave a zero denominator undefined.

Early-year poverty shares are held at their 2008 grid distribution for 2004–2007. Later temporal estimates and allocation methods are described in the parent README. The data cover 2004–2023.
'''
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, data):
    temporary = path.with_suffix(path.suffix+'.partial')
    temporary.write_text(json.dumps(data,indent=2),encoding='utf-8')
    os.replace(temporary,path)


def member(archive, ending):
    matches = [n for n in archive.namelist() if n.replace('\\','/').endswith(ending)]
    require(len(matches)==1, 'Missing or ambiguous ZIP member: '+ending)
    return matches[0]


def load_series(archive, category, stem, prefix):
    name = member(archive, f'/ACS 5-yr estimates/{category}/{stem}.csv')
    with archive.open(name) as stream:
        header = pd.read_csv(stream,nrows=0).columns
    values = [c for c in header if re.fullmatch(prefix+r'[A-Z]{2}\d{2}5',c)]
    keys = [f'GJOIN{year}' for year in range(2010,2024)]
    parts = []
    with archive.open(name) as stream:
        for part in pd.read_csv(stream,usecols=['STATEFP']+keys+values,
                                dtype={k:str for k in ['STATEFP']+keys},chunksize=5000):
            parts.append(part[part.STATEFP.isin(['24','36','39','42','54'])])
    result = pd.concat(parts,ignore_index=True)
    require(len(result)>0,'No relevant states found in '+name)
    return result


def selected(frame, key, columns):
    result = frame[[key]+columns].dropna(subset=[key]).copy()
    result[key] = result[key].astype(str).str.strip()
    require(not result[key].duplicated().any(), 'Duplicate geography IDs: '+key)
    return result.set_index(key)


def row_sum(frame, columns):
    return frame[columns].sum(axis=1,min_count=len(columns))


def age_groups(frame, columns, suffix=None):
    output = pd.DataFrame(index=frame.index)
    if suffix is not None:
        # NHGIS B58 has 18 disjoint age bins for each sex.
        codes = sorted(c for c in columns if c.endswith(suffix))
        require(len(codes)==36, 'Unexpected ACS age schema.')
        for offset, sex in [(0,'male'),(18,'female')]:
            output[sex+'_under5'] = row_sum(frame,codes[offset:offset+1])
            output[sex+'_5_64'] = row_sum(frame,codes[offset+1:offset+15])
            output[sex+'_65plus'] = row_sum(frame,codes[offset+15:offset+18])
    else:
        codes = [f'FMZ{i:03}' for i in range(1,47)]
        for offset, sex in [(0,'male'),(23,'female')]:
            output[sex+'_under5'] = row_sum(frame,codes[offset:offset+1])
            output[sex+'_5_64'] = row_sum(frame,codes[offset+1:offset+17])
            output[sex+'_65plus'] = row_sum(frame,codes[offset+17:offset+23])
    return output[AGE]


def raw_sources(archive):
    age = load_series(archive,'Sex and age','nhgis0018_ts_nominal_tract','B58')
    race = load_series(archive,'Race and ethnicity','nhgis0019_ts_nominal_tract','CV4')
    for end in range(2010,2024):
        suffix = f'{end%100:02}5'
        acols = sorted(c for c in age if re.fullmatch(r'B58[A-Z]{2}'+suffix,c))
        rcols = sorted(c for c in race if re.fullmatch(r'CV4[A-Z]{2}'+suffix,c))
        require(len(acols)==36 and len(rcols)==10,'Unexpected ACS category count.')
        a = selected(age,f'GJOIN{end}',acols)
        r = selected(race,f'GJOIN{end}',rcols)
        result = age_groups(a,acols,suffix)
        # CV4's five non-Hispanic groups and five Hispanic counterparts partition total population.
        for key, cols in zip(RACE,[[rcols[0]],[rcols[1]],[rcols[2]],rcols[3:5],rcols[5:10]]):
            result[key] = row_sum(r,cols).reindex(result.index)
        books = [n for n in archive.namelist() if 'Ratio of Income to Poverty' in n and n.endswith(f'_{end}5_tract_codebook.txt')]
        require(len(books)==1,'Missing income codebook for '+str(end))
        codebook = books[0]
        text = archive.read(codebook).decode('utf-8-sig')
        fields = re.findall(r'^\s*([A-Z0-9]+E00[1-8]):\s*([^\r\n]+)',text,re.M)
        fields = sorted(fields)
        require(len(fields)==8 and [s.strip() for _,s in fields]==
                ['Total','Under .50','.50 to .99','1.00 to 1.24','1.25 to 1.49','1.50 to 1.84','1.85 to 1.99','2.00 and over'],
                'Income definitions differ from expected categories.')
        icols = [c for c,_ in fields]
        with archive.open(codebook.replace('_codebook.txt','.csv')) as stream:
            inc = pd.read_csv(stream,usecols=['GISJOIN']+icols,dtype={'GISJOIN':str})
        inc = selected(inc,'GISJOIN',icols)
        require(LOW_INCOME_RATIO==2.0,'This version implements the confirmed ratio <2.00 definition only.')
        result['low_income'] = row_sum(inc,icols[1:7]).reindex(result.index)
        result['other_poverty_determined'] = inc[icols[7]].reindex(result.index)
        result['poverty_not_determined'] = result[AGE].sum(axis=1,min_count=6)-inc[icols[0]].reindex(result.index)
        result.index.name='GISJOIN'
        yield end-2,end,result[GROUPS]
    acols = [f'FMZ{i:03}' for i in range(1,47)]
    rcols = [f'FMS{i:03}' for i in range(1,15)]
    with archive.open(member(archive,'/Decennial/Sex and age/2000_tract.csv')) as stream:
        a = pd.read_csv(stream,usecols=['GISJOIN']+acols,dtype={'GISJOIN':str})
    with archive.open(member(archive,'/Decennial/Race and ethnicity/2000_tract.csv')) as stream:
        r = pd.read_csv(stream,usecols=['GISJOIN']+rcols,dtype={'GISJOIN':str})
    a,r = selected(a,'GISJOIN',acols),selected(r,'GISJOIN',rcols)
    result = age_groups(a,acols)
    for key, cols in zip(RACE,[[rcols[0]],[rcols[1]],[rcols[2],rcols[3],rcols[4]],rcols[5:7],rcols[7:14]]):
        result[key] = row_sum(r,cols).reindex(result.index)
    # No 2000 poverty-ratio anchor is used. Early income shares use ACS 2006-2010.
    for key in INCOME:
        result[key] = 0.0
    result.index.name='GISJOIN'
    yield 2000,2000,result[GROUPS]


def allocate(task):
    year,end,raw,ids,overlap,expected_hash,cache,cache_hash = task
    require(digest(overlap)==expected_hash,'Overlap file checksum differs: '+str(overlap))
    require(digest(cache)==cache_hash,'Population anchor cache checksum differs: '+str(cache))
    links = pd.read_csv(overlap,dtype={'GISJOIN':str})
    used = pd.Index(links.GISJOIN.unique())
    values = raw.reindex(used).to_numpy(dtype=float)
    require(np.isfinite(values).all() and (values>=0).all(),f'Missing or negative demographic source cells for anchor {year}.')
    total = values[:,:6].sum(axis=1)
    require(np.allclose(values[:,6:11].sum(axis=1),total,rtol=0,atol=1e-6),f'Age and race source totals disagree in {year}.')
    if year!=2000:
        require(np.allclose(values[:,11:14].sum(axis=1),total,rtol=0,atol=1e-6),f'Poverty categories do not reconcile in {year}.')
    positions = np.searchsorted(ids,links.grid_ID.to_numpy())
    require((positions<len(ids)).all() and np.array_equal(ids[positions],links.grid_ID.to_numpy()),'Overlap grid IDs differ.')
    source_index = used.get_indexer(links.GISJOIN)
    weights = links.area_weight.to_numpy()
    require(np.isfinite(weights).all() and (weights>0).all() and (weights<=1+1e-8).all(),'Invalid overlap weights.')
    out = np.zeros((len(ids),len(GROUPS)),dtype=float)
    np.add.at(out,positions,values[source_index]*weights[:,None])
    expected = np.load(cache,allow_pickle=False)
    require(expected.shape==(len(ids),) and np.allclose(out[:,:6].sum(axis=1),expected,rtol=1e-9,atol=1e-6),
            f'Allocated demographics do not reproduce the existing total-population anchor {year}.')
    return year,out


def fractions(values):
    totals = values.sum(axis=1)
    return np.divide(values,totals[:,None],out=np.zeros_like(values),where=totals[:,None]>0)


def extrapolate(recent,year):
    last = recent[-1]
    x = np.arange(-4,1,dtype=float)[:,None,None]
    slope = (x*(recent-last)).sum(axis=0)/30.0
    slope = np.clip(slope,-last/2,last/4)
    return np.maximum(last+(year-2021)*slope,0)


def round_partition(weights,total):
    require(np.isfinite(weights).all() and (weights>=0).all(),'Invalid category shares.')
    sums=weights.sum(axis=1)
    require(((total==0)|(sums>0)).all(),'Positive population has no demographic distribution.')
    values = fractions(weights)*total[:,None]
    result = np.floor(values).astype(np.int64)
    remainder = total-result.sum(axis=1)
    require((remainder>=0).all() and (remainder<=weights.shape[1]).all(),'Invalid partition rounding remainder.')
    order = np.argsort(-(values-result),axis=1,kind='stable')
    for rank in range(weights.shape[1]):
        rows=np.where(remainder>rank)[0]
        result[rows,order[rows,rank]]+=1
    require(np.array_equal(result.sum(axis=1),total) and (result>=0).all(),'Rounded categories do not match total population.')
    return result


def annual_weights(anchors,year):
    if year<2008:
        w=(year-2000)/8
        result=(1-w)*anchors[2000]+w*anchors[2008]
        # Poverty-ratio counts were not supplied for a compatible 2000 anchor.
        result[:,11:14]=fractions(anchors[2008][:,11:14])
        return result
    if year<=2021:
        return anchors[year]
    result=extrapolate(np.stack([anchors[y] for y in range(2017,2022)]),year)
    for start,stop in [(0,6),(6,11),(11,14)]:
        empty=result[:,start:stop].sum(axis=1)==0
        result[empty,start:stop]=anchors[2021][empty,start:stop]
    return result


def output_frames(year,frame,weights):
    ids=frame.grid_ID.to_numpy(dtype=np.int64)
    total=frame.population.to_numpy(dtype=np.int64)
    age=round_partition(weights[:,:6],total)
    race=round_partition(weights[:,6:11],total)
    income_weights=weights[:,11:14].copy()
    unavailable=(total>0)&(income_weights.sum(axis=1)==0)
    require(not unavailable.any() or year<2008,
            f'Unexpected absent income distribution in {year}; only early-year reference gaps can be retained as unavailable.')
    # A zero-population reference supplies no income shares. Preserve these
    # residents explicitly as unclassified, without inventing an income split.
    income_weights[unavailable,2]=1.0
    income=round_partition(income_weights,total)
    coordinates={c:frame[c].to_numpy() for c in ['grid_ID','latitude','longitude']}
    under=age[:,0]+age[:,3]
    over=age[:,2]+age[:,5]
    male=age[:,:3].sum(axis=1);female=age[:,3:].sum(axis=1)
    require((under+over<=total).all() and np.array_equal(male+female,total),'Age or sex partition mismatch.')
    tables={
      f'Age_less_5/Age_less_then_5_{year}.csv':pd.DataFrame(dict(**coordinates,under5_population=under)),
      f'Age_more_64/Age_over_64_{year}.csv':pd.DataFrame(dict(**coordinates,age65_population=over)),
      f'Gender/Sex_{year}.csv':pd.DataFrame(dict(**coordinates,male_population=male,female_population=female)),
      f'Race-ethnicity/Race_{year}.csv':pd.DataFrame(dict(**coordinates,**{k:race[:,i] for i,k in enumerate(RACE)})),
      f'Low_Income/Population_low_income_{year}.csv':pd.DataFrame(dict(**coordinates,low_income_population=income[:,0],
        other_poverty_determined_population=income[:,1],poverty_not_determined_population=income[:,2],
        poverty_determined_population=income[:,0]+income[:,1]))}
    summary=dict(year=year,grid_count=len(ids),total_population=int(total.sum()),under5=int(under.sum()),age65plus=int(over.sum()),
                 male=int(male.sum()),female=int(female.sum()),**{k:int(race[:,i].sum()) for i,k in enumerate(RACE)},
                 low_income=int(income[:,0].sum()),poverty_determined=int(income[:,:2].sum()),
                 poverty_not_determined=int(income[:,2].sum()),
                 income_distribution_unavailable_cells=int(unavailable.sum()),
                 income_distribution_unavailable_people=int(total[unavailable].sum()),
                 all_cell_constraints_passed=True)
    return tables,summary


def main():
    import fcntl
    require(CENSUS_ZIP.is_file(),'Cannot find the existing census ZIP: '+str(CENSUS_ZIP))
    require((GENERAL/'BUILD_SUMMARY.json').is_file(),'Cannot find the installed total-population summary.')
    require((POPULATION_BUILD/'working/checkpoint.sqlite3').is_file(),
            'Keep the total-population working folder; its checkpoint and overlap files are required.')
    OUTPUT_DIR.mkdir(parents=True,exist_ok=True)
    working=OUTPUT_DIR/'working';working.mkdir(exist_ok=True)
    lock=(working/'build.lock').open('a')
    try:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError('Another demographic build is using this output folder.')
    pop_summary=json.loads((GENERAL/'BUILD_SUMMARY.json').read_text())
    require(pop_summary['state']=='complete_total_population','Installed total population is incomplete.')
    expected={r['year']:r for r in pop_summary['annual']}
    require(sorted(expected)==list(range(2004,2024)),'Missing total-population years.')
    first=pd.read_csv(GENERAL/'population_2008.csv')
    ids=first.grid_ID.to_numpy(dtype=np.int64)
    require(len(ids)==pop_summary['grid_cells'] and (np.diff(ids)>0).all(),'Unexpected grid IDs.')
    for year in range(2004,2024):
        require(digest(GENERAL/f'population_{year}.csv')==expected[year]['file_sha256'],'Installed population file checksum differs.')
    with zipfile.ZipFile(CENSUS_ZIP) as archive:
        zip_identity=hashlib.sha256(json.dumps([(i.filename,i.CRC,i.file_size) for i in archive.infolist()]).encode()).hexdigest()
    identity=dict(version=VERSION,population_signature=pop_summary['signature'],census_zip_identity=zip_identity,
                  population_file_hashes={y:expected[y]['file_sha256'] for y in sorted(expected)},low_income_ratio=LOW_INCOME_RATIO,
                  race_groups=RACE,early_income='Hold 2008 grid shares for 2004-2007',rounding='Joint age/sex and separate race and poverty partitions by grid')
    signature=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    identity_path=working/'input_identity.json'
    if identity_path.exists():
        require(json.loads(identity_path.read_text())['signature']==signature,'Inputs changed; select a new OUTPUT_DIR.')
    atomic_json(identity_path,dict(signature=signature,**identity))
    db=sqlite3.connect(working/'checkpoint.sqlite3')
    db.execute('CREATE TABLE IF NOT EXISTS products (name TEXT PRIMARY KEY,sha256 TEXT)')
    base_db=sqlite3.connect(POPULATION_BUILD/'working/checkpoint.sqlite3')
    anchor_hashes=dict(base_db.execute('SELECT name,sha256 FROM products'));base_db.close()
    anchors={};pending=[]
    with zipfile.ZipFile(CENSUS_ZIP) as archive:
        for year,end,raw in raw_sources(archive):
            name=f'anchor_{year}.npy';path=working/name
            old=db.execute('SELECT sha256 FROM products WHERE name=?',(name,)).fetchone()
            if path.exists():
                require(old and digest(path)==old[0],'Unverified demographic anchor cache.')
                anchors[year]=np.load(path,allow_pickle=False)
                print(f'Reused demographic anchor {year}',flush=True)
                continue
            info=json.loads((POPULATION_BUILD/f'working/overlaps/summary_{end}.json').read_text())
            cache=POPULATION_BUILD/f'working/anchor_{year}.npy'
            require(cache.name in anchor_hashes,'Missing total-population anchor checkpoint.')
            pending.append((year,end,raw,ids,POPULATION_BUILD/'working/overlaps'/info['output_file'],info['output_sha256'],cache,anchor_hashes[cache.name]))
    workers=min(WORKERS,int(os.environ.get('SLURM_CPUS_PER_TASK',WORKERS)))
    require(workers>0,'WORKERS must be positive.')
    print(f'Allocating {len(pending)} source anchors; workers={workers}',flush=True)
    def save_anchor(result):
        year,values=result;path=working/f'anchor_{year}.npy';tmp=path.with_suffix('.partial')
        with tmp.open('wb') as stream:np.save(stream,values,allow_pickle=False)
        db.execute('INSERT OR REPLACE INTO products VALUES (?,?)',(path.name,digest(tmp)));db.commit();os.replace(tmp,path)
        anchors[year]=values
        print(f'Completed demographic anchor {year}',flush=True)
    if workers==1:
        for task in pending:save_anchor(allocate(task))
    elif pending:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            futures=[pool.submit(allocate,task) for task in pending]
            for future in concurrent.futures.as_completed(futures):save_anchor(future.result())
    require(sorted(anchors)==[2000]+list(range(2008,2022)),'Missing demographic anchor.')
    summaries=[];products={};unavailable_rows=[]
    for year in range(2004,2024):
        frame=pd.read_csv(GENERAL/f'population_{year}.csv')
        require(np.array_equal(frame.grid_ID.to_numpy(),ids),'Population grids differ between years.')
        require(np.isfinite(frame.population).all() and (frame.population>=0).all() and np.equal(frame.population,np.floor(frame.population)).all(),'Invalid total population.')
        weights=annual_weights(anchors,year)
        tables,summary=output_frames(year,frame,weights)
        unavailable=(frame.population.to_numpy()>0)&(weights[:,11:14].sum(axis=1)==0)
        for row in frame.loc[unavailable,['grid_ID','population']].itertuples(index=False):
            unavailable_rows.append(dict(year=year,grid_ID=int(row.grid_ID),population=int(row.population),
                                         reason='2008 reference population is zero; early-year income distribution unavailable'))
        for relative,table in tables.items():
            path=OUTPUT_DIR/'demographics'/relative;path.parent.mkdir(parents=True,exist_ok=True)
            old=db.execute('SELECT sha256 FROM products WHERE name=?',(relative,)).fetchone()
            if path.exists():
                require(old and digest(path)==old[0],'Unverified existing demographic output: '+str(path))
            else:
                tmp=path.with_suffix('.csv.partial');table.to_csv(tmp,index=False,float_format='%.15g')
                db.execute('INSERT OR REPLACE INTO products VALUES (?,?)',(relative,digest(tmp)));db.commit();os.replace(tmp,path)
            products[relative]=digest(path)
        summaries.append(summary)
        print(f'{year}: five demographic files saved/verified; all grid constraints passed',flush=True)
    pd.DataFrame(summaries).to_csv(OUTPUT_DIR/'annual_demographic_summary.csv',index=False)
    availability=OUTPUT_DIR/'demographics/income_distribution_unavailable.csv'
    temporary=availability.with_suffix('.csv.partial')
    pd.DataFrame(unavailable_rows,columns=['year','grid_ID','population','reason']).to_csv(temporary,index=False)
    os.replace(temporary,availability)
    for relative,text in READMES.items():
        path=OUTPUT_DIR/'demographics'/relative
        require(not path.exists() or path.read_text(encoding='utf-8')==text
                or hashlib.sha256(path.read_text(encoding='utf-8').encode('utf-8')).hexdigest()==PREVIOUS_README_HASHES.get(relative),
                'An unexpected README exists and was retained: '+str(path))
        path.parent.mkdir(parents=True,exist_ok=True)
        temporary=path.with_suffix('.md.partial');temporary.write_text(text,encoding='utf-8');os.replace(temporary,path)
    atomic_json(OUTPUT_DIR/'BUILD_SUMMARY.json',dict(state='complete_demographics',signature=signature,
        years=20,grid_cells=len(ids),files=len(products),all_cell_constraints_passed=True,definitions=identity,
        output_policy=OUTPUT_POLICY,income_distribution_unavailable_records=len(unavailable_rows),
        availability_file_sha256=digest(availability),
        limitations=['Early income shares are held at the ACS 2006-2010 grid distribution for 2004-2007; where that reference population is zero, early residents retain unavailable income status.',
          'Recent demographic category trends are extrapolated with nonnegative slope constraints and reconciled to each installed grid population.',
          'The combined non-Hispanic AIAN/Asian/Pacific group cannot identify Asian populations separately.',
          'Demographic category totals are reconciled to integer grid populations; they are not independently controlled official annual totals.',
          'The total-population midpoint convention, early interpolation and recent geometry assumptions also apply.'],
        annual=summaries,file_sha256=products))
    db.close();lock.close()
    print(f'Complete: {OUTPUT_DIR}. Share BUILD_SUMMARY.json, annual_demographic_summary.csv and demographics/income_distribution_unavailable.csv before installing shared files.',flush=True)


if __name__=='__main__':
    main()
