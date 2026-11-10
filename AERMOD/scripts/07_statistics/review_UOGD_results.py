"""Independent scratch statistics for concentrations, contributions or exposure.
Run with one of the three batch files. Existing shared outputs are read only.
"""
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor,as_completed
import argparse,calendar,hashlib,json,os,sqlite3,contextlib,multiprocessing
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np
import pandas as pd
import xarray as xr
ROOT=Path(os.environ.get('TRIMAQS_DATA_ROOT', str(Path(__file__).resolve().parents[3]/'TriMAQs_data')))
PROCESSED=ROOT/'working/statistics_inputs'
CONCENTRATIONS=PROCESSED/'concentrations'
CONTRIBUTIONS=PROCESSED/'source_apportionment'
YEARS=list(range(2004,2024))
BLOCK_DAYS=7
MODE=None
SOURCES=['Tier_1','Tier_2','UOGD']
POLS={'NO2':'NO2','PM25':'PM2.5','VOC':'VOC','SO2':'SO2'}
UNITS={'NO2':'ppb','SO2':'ppb','PM25':'ug/m3','VOC':'ug/m3'}
# Same 2021 guideline values used in the report; compare in the NetCDF units.
# These source-only threshold counts are descriptive, not total-ambient compliance.
ANNUAL={'NO2':10*24.45/46.01,'PM25':5.}
DAILY={'NO2':25*24.45/46.01,'PM25':15.,'SO2':40*24.45/64.07}
POP_SPECS={
 'General':('general/population_{year}.csv',['Population','population','pop_general']),
 'Under5':('demographics/Age_less_5/Age_less_then_5_{year}.csv',['under5_population','Under5_population']),
 'Over64':('demographics/Age_more_64/Age_over_64_{year}.csv',['age65_population','Age65_population']),
 'LowIncome':('demographics/Low_Income/Population_low_income_{year}.csv',['low_income_population','Low_income_population']),
 'White_non_Hispanic':('demographics/Race-ethnicity/Race_{year}.csv',['white_non_hispanic']),
 'Black_non_Hispanic':('demographics/Race-ethnicity/Race_{year}.csv',['black_non_hispanic']),
 'AIAN_Asian_Pacific_non_Hispanic':('demographics/Race-ethnicity/Race_{year}.csv',['aian_asian_pacific_non_hispanic']),
 'Hispanic':('demographics/Race-ethnicity/Race_{year}.csv',['hispanic']),
 'Other_non_Hispanic':('demographics/Race-ethnicity/Race_{year}.csv',['other_non_hispanic']),
 'Male':('demographics/Gender/Sex_{year}.csv',['male_population']),
 'Female':('demographics/Gender/Sex_{year}.csv',['female_population'])}
GROUP_LABELS={
 'General':'Total population','Under5':'Age under 5','Over64':'Age 65 and over',
 'LowIncome':'Income below twice the poverty threshold',
 'White_non_Hispanic':'Non-Hispanic White','Black_non_Hispanic':'Non-Hispanic Black',
 'AIAN_Asian_Pacific_non_Hispanic':'Non-Hispanic American Indian/Alaska Native, Asian and Pacific Islander (combined)',
 'Hispanic':'Hispanic or Latino (any race)','Other_non_Hispanic':'Other non-Hispanic races and multiracial',
 'Male':'Male','Female':'Female'}


def require(ok,message):
    if not ok:raise ValueError(message)

def frame_values(frame,ids,candidates):
    lower={str(c).lower():c for c in frame.columns}
    require('grid_id' in lower,'Population grid_ID column missing')
    key=lower['grid_id'];frame=frame.copy();frame[key]=pd.to_numeric(frame[key],errors='raise').astype('int64')
    require(not frame[key].duplicated().any(),'Duplicate population grid IDs')
    column=next((lower[c.lower()] for c in candidates if c.lower() in lower),None)
    require(column is not None,'Population count field missing: '+str(candidates))
    values=pd.to_numeric(frame.set_index(key)[column],errors='raise').reindex(ids).to_numpy(dtype=float)
    require(np.isfinite(values).all() and (values>=0).all(),'Missing/invalid population cells; no zeros have been substituted')
    return values

def population(year,ids):
    result={};issues=[]
    for group,(pattern,columns) in POP_SPECS.items():
        file=ROOT/'inputs/population_demographics'/pattern.format(year=year)
        if not file.is_file():issues.append(dict(year=year,group=group,issue='not supplied',file=str(file)));continue
        try:
            values=frame_values(pd.read_csv(file),ids,columns)
            require(group!='General' or values.sum()>0,'General population sum is zero')
            result[group]=values
        except Exception as e:issues.append(dict(year=year,group=group,issue=str(e),file=str(file)))
    return result,issues

def distribution(values):
    values=np.asarray(values,dtype=float);valid=values[np.isfinite(values)]
    return dict(valid_cells=len(valid),mean=float(valid.mean()) if len(valid) else None,
                median=float(np.median(valid)) if len(valid) else None,
                p95=float(np.quantile(valid,.95)) if len(valid) else None,
                p99=float(np.quantile(valid,.99)) if len(valid) else None,
                maximum=float(valid.max()) if len(valid) else None)

def rrr_rows(year,pol,kind,pop,amount,days=1):
    if 'General' not in pop:return []
    total=float(pop['General'].sum())*days;exposed=float(amount['General'])
    rows=[]
    for group,weights in pop.items():
        denom=float(weights.sum())*days
        ratio=(amount[group]/denom)/(exposed/total) if denom>0 and total>0 and exposed>0 else np.nan
        rows.append(dict(Year=year,Pollutant=POLS[pol],Group=group,Group_label=GROUP_LABELS[group],kind=kind,
                         Total_population=float(weights.sum()),denominator_person_days=denom,
                         Exposed_population=float(amount[group]) if kind=='annual' else None,
                         exposed_person_days=float(amount[group]) if kind=='daily_aggregated' else None,RRR=float(ratio)))
    return rows

def one_year(year):
    result={k:[] for k in ['annual','daily','seasons','fractions','exposure','annual_thresholds','daily_thresholds','rrr','population_issues','old_comparison','old_fraction_comparison','income_availability']}
    with xr.open_dataset(CONCENTRATIONS/f'aermod_daily_conc_{year}.nc') as ds:
        require(list(ds.source.values)==SOURCES,'Unexpected source labels')
        mask=ds.study_mask.values.astype(bool);inside=~mask;ids=ds.grid_ID.values[inside].astype('int64')
        lat=ds.lat.values[inside];lon=ds.lon.values[inside]
        dates=pd.DatetimeIndex(ds.time.values)
        require(dates.equals(pd.date_range(f'{year}-01-01',f'{year}-12-31')),'Incomplete calendar')
        pop={}
        if MODE=='exposure':
            pop,issues=population(year,ids);result['population_issues']=issues
            require(not issues,'Population inputs are incomplete: '+str(issues[:3]))
            availability_file=ROOT/'inputs/population_demographics/demographics/income_distribution_unavailable.csv'
            require(availability_file.is_file(),'Missing demographic income-availability report.')
            availability=pd.read_csv(availability_file)
            selected=availability[availability.year==year]
            require(not selected.duplicated('grid_ID').any() and selected.grid_ID.isin(ids).all(),'Invalid income-availability grid IDs.')
            result['income_availability'].append(dict(year=year,unavailable_cells=len(selected),
                unavailable_people=int(selected.population.sum()),note='Income status unavailable in zero-population reference cells; residents retained in total population.'))
        if MODE in {'concentrations','exposure'}:
            for pol in POLS:
                require(ds[pol].attrs.get('units')==UNITS[pol],'Units differ for '+pol)
                annual=ds[pol+'_annual'].transpose('source','y','x').values[:,inside].astype(float)
                require(np.isfinite(annual).all(),'Missing annual concentrations inside study boundary')
                means=[];daily_exposed={g:0. for g in pop};maximum=np.full(3,-np.inf);maxdate=['']*3;maxgrid=[None]*3
                for start in range(0,len(dates),BLOCK_DAYS):
                    values=ds[pol].isel(time=slice(start,start+BLOCK_DAYS)).transpose('time','source','y','x').values[:,:,inside].astype(float)
                    require(np.isfinite(values).all() and (values>=0).all(),'Missing/negative daily concentration inside study boundary')
                    require(np.allclose(values[:,0]+values[:,1],values[:,2],rtol=2e-5,atol=1e-8),'Tier sum differs')
                    blockmean=values.mean(axis=2);means.append(blockmean)
                    for si,source in enumerate(SOURCES):
                        block=values[:,si,:];at=np.unravel_index(np.argmax(block),block.shape)
                        if block[at]>maximum[si]:maximum[si]=block[at];maxdate[si]=str(dates[start+at[0]].date());maxgrid[si]=int(ids[at[1]])
                        result['daily'] += [dict(year=year,date=str(date.date()),pollutant=POLS[pol],source=source,units=UNITS[pol],spatial_mean=float(v)) for date,v in zip(dates[start:start+len(blockmean)],blockmean[:,si])]
                    if pol in DAILY:
                        above=values[:,2,:]>DAILY[pol]
                        for j,date in enumerate(dates[start:start+len(values)]):
                            result['daily_thresholds'].append(dict(Date=str(date.date()),Pollutant=POLS[pol],threshold=DAILY[pol],units=UNITS[pol],grids_above=int(above[j].sum()),Population_Above=float(pop['General'][above[j]].sum()) if 'General' in pop else None))
                        counts=above.sum(axis=0)
                        for group,weights in pop.items():daily_exposed[group]+=float(np.dot(counts,weights))
                dailymeans=np.concatenate(means)
                for si,source in enumerate(SOURCES):
                    a=annual[si];at=int(np.argmax(a));stats=distribution(a)
                    require(np.isclose(a.mean(),dailymeans[:,si].mean(),rtol=2e-5,atol=1e-8),'Annual mean differs from daily mean')
                    result['annual'].append(dict(year=year,pollutant=POLS[pol],source=source,units=UNITS[pol],**stats,
                        annual_max_grid_id=int(ids[at]),annual_max_lat=float(lat[at]),annual_max_lon=float(lon[at]),
                        maximum_daily_grid_concentration=float(maximum[si]),daily_grid_max_date=maxdate[si],daily_grid_max_id=maxgrid[si],
                        maximum_daily_spatial_mean=float(dailymeans[:,si].max()),daily_spatial_mean_peak_date=str(dates[int(np.argmax(dailymeans[:,si]))].date())))
                    for season,months in {'winter':[12,1,2],'spring':[3,4,5],'summer':[6,7,8],'fall':[9,10,11]}.items():
                        v=float(dailymeans[np.isin(dates.month,months),si].mean());base=float(a.mean())
                        result['seasons'].append(dict(year=year,pollutant=POLS[pol],source=source,season=season,seasonal_mean=v,annual_mean=base,seasonal_factor=v/base if base>0 else None,definition='calendar-year months; December is in the same calendar year'))
                try:result['old_comparison'] += compare_previous(year,pol,ids,annual[2])
                except Exception as e:result['old_comparison'].append(dict(year=year,pollutant=POLS[pol],state='old_baseline_needs_review',error=str(e)))
                for group,weights in pop.items():
                    denominator=weights.sum();v=float(np.dot(annual[2],weights)/denominator) if denominator>0 else np.nan
                    result['exposure'].append(dict(Year=year,Pollutant=POLS[pol],Group=group,Group_label=GROUP_LABELS[group],units=UNITS[pol],population=float(denominator),annual_spatial_mean=float(annual[2].mean()),weighted_concentration=v,exposure_concentration_ratio=v/annual[2].mean() if annual[2].mean()>0 else None))
                if pol in ANNUAL:
                    above=annual[2]>ANNUAL[pol];amount={g:float(weights[above].sum()) for g,weights in pop.items()}
                    result['annual_thresholds'].append(dict(Year=year,Pollutant=POLS[pol],threshold=ANNUAL[pol],units=UNITS[pol],grids_above=int(above.sum()),Population_Above=amount.get('General')))
                    result['rrr']+=rrr_rows(year,pol,'annual',pop,amount)
                if pol in DAILY:result['rrr']+=rrr_rows(year,pol,'daily_aggregated',pop,daily_exposed,len(dates))
    if MODE=='source_contributions':
        with xr.open_dataset(CONTRIBUTIONS/f'uogd_source_apportionment_{year}.nc') as ds:
            require(list(ds.source.values)==SOURCES and np.array_equal(ds.study_mask.values.astype(bool),mask),'Source fraction grid differs')
            require(np.array_equal(ds.grid_ID.values[inside],ids) and np.allclose(ds.lat.values[inside],lat) and np.allclose(ds.lon.values[inside],lon),'Source fraction coordinates differ')
            for pol in ['NO2','PM25','SO2']:
                days=ds[pol+'_valid_days'].values[inside]
                for method in ['annual_mean_fraction','annual_ratio']:
                    require(ds[pol+'_'+method].attrs.get('units')=='1','Source-ratio units differ')
                    arrays=ds[pol+'_'+method].transpose('source','y','x').values[:,inside]
                    require(np.array_equal(np.isfinite(arrays),np.broadcast_to(days>0,arrays.shape)) and (arrays[np.isfinite(arrays)]>=0).all(),'Annual ratio/valid-day mask differs')
                    result['old_fraction_comparison'] += compare_previous(year,pol,ids,arrays[2],method)
                    for si,source in enumerate(SOURCES):
                        v=arrays[si];result['fractions'].append(dict(year=year,pollutant=POLS[pol],source=source,statistic=method,**distribution(v),
                            cells_at_or_above_one=int((v>=1).sum()),cells_above_003=int((v>.03).sum()),cells_above_008=int((v>.08).sum()),cells_above_040=int((v>.4).sum()),
                            study_cells_no_valid_days=int((days==0).sum()),minimum_valid_days=int(days.min()),median_valid_days=float(np.median(days)),maximum_valid_days=int(days.max())))
    return result


def compare_previous(year,pol,ids,current,method=None):
    fraction=method is not None
    path=ROOT/(f'outputs/source_contributions/uogd_source_apportionment_{year}.nc' if fraction else f'outputs/concentrations/aermod_daily_conc_{year}.nc')
    result=dict(year=year,pollutant=POLS[pol],baseline_file=str(path),new_statistic=method or 'annual_mean_concentration')
    if not path.is_file():return [dict(result,state='previous_shared_file_absent')]
    try:
        with xr.open_dataset(path) as old:
            require(list(old.source.values)==SOURCES,'Previous source labels differ')
            inside=~old.study_mask.values.astype(bool);oldids=old.grid_ID.values[inside].astype(np.int64)
            require(len(np.unique(oldids))==len(oldids),'Previous file has duplicate grid IDs')
            name=pol+'_'+method if fraction else pol+'_annual'
            if not fraction:require(old[name].attrs.get('units')==UNITS[pol],'Previous units differ')
            previous=pd.Series(old[name].transpose('source','y','x').values[2,inside],index=oldids).reindex(ids).to_numpy(float)
            oldlat=pd.Series(old.lat.values[inside],index=oldids).reindex(ids).to_numpy(float)
            oldlon=pd.Series(old.lon.values[inside],index=oldids).reindex(ids).to_numpy(float)
        with xr.open_dataset(CONCENTRATIONS/f'aermod_daily_conc_{year}.nc') as current_file:
            current_mask=~current_file.study_mask.values.astype(bool)
            current_ids=current_file.grid_ID.values[current_mask].astype(np.int64)
            lat=pd.Series(current_file.lat.values[current_mask],index=current_ids).reindex(ids).to_numpy(float)
            lon=pd.Series(current_file.lon.values[current_mask],index=current_ids).reindex(ids).to_numpy(float)
        common=np.isfinite(oldlat)&np.isfinite(oldlon)
        require(np.allclose(oldlat[common],lat[common],rtol=0,atol=1e-5) and np.allclose(oldlon[common],lon[common],rtol=0,atol=1e-5),'Previous grid IDs do not align geographically')
        valid=np.isfinite(previous)&np.isfinite(current)
        if not valid.any():return [dict(result,state='no_common_valid_cells')]
        a=previous[valid];b=current[valid]
        return [dict(result,state='compared',common_valid_grids=int(valid.sum()),units='1' if fraction else UNITS[pol],
            old_mean_common_cells=float(a.mean()),new_mean_common_cells=float(b.mean()),mean_difference=float((b-a).mean()),
            percent_change_of_mean=float(100*(b.mean()/a.mean()-1)) if a.mean()!=0 else None,
            old_maximum=float(a.max()),new_maximum=float(b.max()),grid_rmse=float(np.sqrt(np.mean((b-a)**2))))]
    except Exception as exc:return [dict(result,state='previous_baseline_needs_review',error=str(exc))]


TABLES={
    'concentrations':['annual','daily','seasons','old_comparison'],
    'source_contributions':['fractions','old_fraction_comparison'],
    'exposure':['exposure','annual_thresholds','daily_thresholds','rrr','population_issues','income_availability'],
}


def year_fingerprint(year):
    files=[CONCENTRATIONS/f'aermod_daily_conc_{year}.nc']
    if MODE=='source_contributions':files += [CONTRIBUTIONS/f'uogd_source_apportionment_{year}.nc',ROOT/f'outputs/source_contributions/uogd_source_apportionment_{year}.nc']
    if MODE=='concentrations':files += [ROOT/f'outputs/concentrations/aermod_daily_conc_{year}.nc']
    if MODE=='exposure':
        files += sorted(set(ROOT/'inputs/population_demographics'/p.format(year=year) for p,_ in POP_SPECS.values()))
        files += [ROOT/'inputs/population_demographics/demographics/income_distribution_unavailable.csv']
    payload=dict(mode=MODE,year=year,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),files=[])
    for path in files:
        item=dict(path=str(path.resolve()),exists=path.is_file())
        if path.is_file():s=path.stat();item.update(size=s.st_size,mtime_ns=s.st_mtime_ns)
        payload['files'].append(item)
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()


def setup(config):
    global ROOT,PROCESSED,CONCENTRATIONS,CONTRIBUTIONS,MODE,BLOCK_DAYS
    ROOT=Path(config['root']);PROCESSED=Path(config['processed']);MODE=config['mode'];BLOCK_DAYS=config['block_days']
    CONCENTRATIONS=PROCESSED/'concentrations';CONTRIBUTIONS=PROCESSED/'source_apportionment'


def task(year):
    return one_year(year)


def write_stage(results,destination):
    for name in TABLES[MODE]:
        records=results.get(name,[])
        fields=['year','group','issue','file'] if name=='population_issues' else None
        frame=pd.DataFrame(records,columns=fields) if not records and fields else pd.DataFrame(records)
        partial=destination/(name+'.partial.csv');frame.to_csv(partial,index=False);partial.replace(destination/(name+'.csv'))
    if MODE=='exposure':
        pd.DataFrame(results['exposure']).to_csv(destination/'annual_exposure.csv',index=False)
        pd.DataFrame(results['annual_thresholds']).pivot(index='Year',columns='Pollutant',values='Population_Above').reset_index().to_csv(destination/'Population_Above_Limit_Annual.csv',index=False)
        pd.DataFrame(results['daily_thresholds']).to_csv(destination/'Population_Above_Limit_Daily.csv',index=False)
        ratios=pd.DataFrame(results['rrr'])
        ratios[ratios.kind=='annual'].to_csv(destination/'RRR_annual_allPollutants.csv',index=False)
        for pol in ['NO2','PM2.5','SO2']:ratios[(ratios.kind=='daily_aggregated')&(ratios.Pollutant==pol)].to_csv(destination/f'RRR_dailyAnnual_{pol}.csv',index=False)
        pd.DataFrame([dict(Group=k,Group_label=v) for k,v in GROUP_LABELS.items()]).to_csv(destination/'population_group_labels.csv',index=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=list(TABLES),required=True)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--processed',type=Path,default=PROCESSED)
    parser.add_argument('--years',type=int,nargs='+',default=YEARS)
    parser.add_argument('--block-days',type=int,default=7)
    args=parser.parse_args()
    require(1<=args.workers<=int(os.environ.get('SLURM_CPUS_PER_TASK','1')),'Workers exceed allocated CPUs')
    require(args.block_days>=1 and len(set(args.years))==len(args.years),'Invalid block size or duplicate years')
    config=dict(root=str(args.root.resolve()),processed=str(args.processed.resolve()),mode=args.stage,block_days=args.block_days)
    setup(config)
    output=PROCESSED/'report_review'/MODE;output.mkdir(parents=True,exist_ok=True)
    import fcntl
    lock=(output/'statistics.lock').open('a+')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise RuntimeError('This statistics stage is already running')
    stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f');logpath=output/('STATISTICS_'+stamp+'.jsonl')
    summarypath=output/('STATISTICS_'+stamp+'_SUMMARY.json');failures=[];done=0
    with contextlib.closing(sqlite3.connect(output/'checkpoint.sqlite3')) as db,logpath.open('x') as log:
        db.execute('CREATE TABLE IF NOT EXISTS years(year INTEGER PRIMARY KEY,fingerprint TEXT,result TEXT)')
        pending=[]
        for year in args.years:
            require((CONCENTRATIONS/f'aermod_daily_conc_{year}.nc').is_file(),'New concentration file absent: '+str(year))
            if MODE=='source_contributions':require((CONTRIBUTIONS/f'uogd_source_apportionment_{year}.nc').is_file(),'New contribution file absent: '+str(year))
            sig=year_fingerprint(year);saved=db.execute('SELECT fingerprint FROM years WHERE year=?',(year,)).fetchone()
            if saved and saved[0]==sig:
                done+=1;log.write(json.dumps(dict(year=year,state='skipped_verified_complete'))+'\n')
            else:pending.append((year,sig))
        with ProcessPoolExecutor(max_workers=args.workers,initializer=setup,initargs=(config,),mp_context=multiprocessing.get_context('spawn')) as pool:
            jobs={pool.submit(task,year):(year,sig) for year,sig in pending}
            for future in as_completed(jobs):
                year,sig=jobs[future]
                try:
                    result=future.result();require(year_fingerprint(year)==sig,'Inputs changed during calculation')
                    with db:db.execute('INSERT OR REPLACE INTO years VALUES (?,?,?)',(year,sig,json.dumps(result,allow_nan=True)))
                    row=dict(year=year,state='completed')
                except Exception as exc:
                    row=dict(year=year,state='failed',error=str(exc));failures.append(row)
                done+=1;log.write(json.dumps(row)+'\n');log.flush();print(year,row['state'],row.get('error',''),flush=True)
                temporary=summarypath.with_suffix('.tmp')
                temporary.write_text(json.dumps(dict(stage=MODE,selected_years=args.years,finished_years=done,failed_years=len(failures),all_complete=False,output=str(output),log=str(logpath)),indent=2));temporary.replace(summarypath)
        require(not failures,'Some years failed; completed checkpoints retained. Share '+str(logpath))
        results={key:[] for key in TABLES[MODE]}
        for year in args.years:
            saved=db.execute('SELECT fingerprint,result FROM years WHERE year=?',(year,)).fetchone()
            require(saved and saved[0]==year_fingerprint(year),'Year checkpoint absent/stale')
            result=json.loads(saved[1])
            for key in results:results[key].extend(result.get(key,[]))
        write_stage(results,output)
    summary=dict(stage=MODE,selected_years=args.years,finished_years=done,failed_years=0,all_complete=True,output=str(output),log=str(logpath))
    summarypath.write_text(json.dumps(summary,indent=2));print('Share CSVs and:',summarypath,flush=True)

if __name__=='__main__':main()
