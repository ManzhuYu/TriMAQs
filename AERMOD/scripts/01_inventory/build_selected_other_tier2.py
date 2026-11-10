"""Build the selected other Tier 2 monthly inventory from reviewed NEI annual inputs."""
from pathlib import Path
import os,calendar
import numpy as np
import pandas as pd
ROOT=Path(os.environ.get('TRIMAQS_DATA_ROOT',str(Path(__file__).resolve().parents[3]/'TriMAQs_data')))
YEARS=np.arange(2004,2024)
ANCHORS=np.array([2011,2014,2017,2020])
MONTHS=[calendar.month_abbr[m].lower()+'_tonnes' for m in range(1,13)]
def require(ok, message):
    if not ok:
        raise ValueError(message)

def annual_from_pool(pool, production):
    require(np.all(pool[ANCHORS-2004] > 0) and np.all(production[ANCHORS-2004] > 0), 'Missing/zero NEI anchor or production value.')
    logef = np.log(pool[ANCHORS-2004] / production[ANCHORS-2004])
    result = []
    for year in YEARS:
        index = max(0, min(len(ANCHORS)-2, int(np.searchsorted(ANCHORS, year))-1))
        t = (year-ANCHORS[index]) / (ANCHORS[index+1]-ANCHORS[index])
        result.append(np.exp(logef[index] + t*(logef[index+1]-logef[index])))
    factors = np.array(result)
    return factors * production, factors

def monthly_from_annual(annual):
    """Reproduce the final archived IDW/Gaussian redistribution, sigma=2.5."""
    values = np.zeros((20, 12))
    midyears = YEARS + 0.5
    gfull = np.exp(-((np.arange(1,13)-6.5)**2)/(2*2.5**2))
    gfull[[0,-1]] = 0
    for i, (year, total) in enumerate(zip(YEARS, annual)):
        if total == 0:
            continue
        for j, time in enumerate(year+(np.arange(1,13)-0.5)/12):
            indices = np.argsort(np.abs(midyears-time))[:2]
            weights = 1/np.abs(time-midyears[indices])
            values[i,j] = np.sum(weights*annual[indices]/12)/weights.sum()
        for _ in range(20):
            diff = total-values[i].sum()
            if abs(diff) < 1e-9:
                break
            weights = gfull.copy()
            weights[values[i] <= 0] = 0
            if weights.sum() == 0:
                weights = (values[i] > 0).astype(float)
            require(weights.sum() > 0, 'Cannot redistribute positive annual emissions.')
            values[i] += weights/weights.sum()*diff
            negative = values[i] < 0
            if negative.any():
                values[i,negative] = 0
            else:
                break
        require(np.isclose(values[i].sum(),total,rtol=1e-10,atol=1e-8), 'Monthly mass balance failed.')
    return values

def main():
 folder=ROOT/'inputs/source_data/other_tier2'
 selection=pd.read_csv(folder/'facility_selection.csv',dtype={'source_id':str})
 raw=pd.read_csv(folder/'nei_annual.csv',dtype={'source_id':str})
 driver=pd.read_csv(folder/'production_driver.csv')
 driver=driver[(driver.kind=='regional_driver')&(driver.state=='read')]
 require(driver.year.is_unique,'Duplicate production year')
 prod=driver.set_index('year').annual_total_as_stored.reindex(YEARS,fill_value=0).to_numpy()
 keep=set(selection.loc[selection.include_in_other_tier2.astype(str).str.lower().eq('true'),'source_id'])
 require(not raw.duplicated(['source_id','year','pollutant']).any(),'Duplicate NEI records')
 output=ROOT/'working/inventory/other_tier2';output.mkdir(parents=True,exist_ok=True)
 for pol in ['NOX','PM2.5','VOC','SO2']:
  selected=raw[raw.source_id.isin(keep)&raw.pollutant.eq(pol)]
  totals=selected.groupby('year').raw_annual_emission.sum().reindex(YEARS,fill_value=0).to_numpy()
  annual,_=annual_from_pool(totals,prod);monthly=monthly_from_annual(annual)
  rows=[]
  for i,y in enumerate(YEARS):
   q=selected[selected.year.eq(y)].copy()
   if annual[i]<=0:continue
   masses=q.raw_annual_emission.to_numpy()[:,None]*(monthly[i]/totals[i] if totals[i]>0 else np.zeros(12))
   require(np.allclose(masses.sum(axis=0),monthly[i],rtol=1e-10,atol=1e-8),'Monthly mass balance failed')
   q=q[['source_id','year','pollutant']].reset_index(drop=True);q[MONTHS]=masses;rows.append(q)
  pd.concat(rows,ignore_index=True).to_csv(output/f'{pol}_monthly.csv',index=False)
 print('Monthly inventory:',output)
if __name__=='__main__':main()
