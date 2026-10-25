#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit

def simulate(n,seed,out,inject=True):
    if n<20:raise ValueError('At least 20 workers required')
    out=Path(out);out.mkdir(exist_ok=True)
    rng=np.random.default_rng(seed);j=np.arange(n)%12;rng.shuffle(j)
    industries=np.array(['Construction','Brick manufacture','Stone fabrication'])
    group=j%3;age=rng.uniform(25,65,n);duration=np.minimum(rng.uniform(3,30,n),age-18)
    site=rng.normal(0,.18,12)
    concentration=np.exp(np.log(np.array([.04,.07,.12])[group])+site[j]+rng.normal(0,.35,n))
    cumulative=concentration*np.maximum(duration-5,0)
    x=np.log2(1+cumulative);smoke=rng.binomial(1,.35,n)
    prevalent=rng.binomial(1,expit(-3.4+.5*x+.3*(age-45)/10),n)
    eta=np.log(1.7)*x+np.log(1.3)*(age-45)/10+np.log(1.5)*smoke
    event_time=(-np.log(rng.uniform(size=n))/(.012*np.exp(eta)))**(1/1.4)
    censor=np.minimum(8,rng.exponential(25,n));time=np.minimum(event_time,censor);event=(event_time<=censor).astype(int)
    workers=pd.DataFrame(dict(worker_id=[f'W{i+1:04d}' for i in range(n)],workplace=[f'S{k+1:02d}' for k in j],industry=industries[group],age=age,smoker=smoke.astype(float),employment_years=duration,historical_mg_m3=concentration,baseline_disease=prevalent))
    follow=pd.DataFrame(dict(worker_id=workers.worker_id,time_years=time,event=event));follow=follow.loc[prevalent==0].copy()
    measurements=pd.DataFrame([dict(worker_id=workers.worker_id[i],sample_id=f'M{i+1:04d}_{k+1}',concentration_mg_m3=concentration[i]*np.exp(rng.normal(0,.25)),duration_hours=8,unit='mg/m3') for i in range(n) for k in range(3)])
    defects=[]
    if inject:
        missing=rng.choice(n,max(1,int(.03*n)),replace=False);workers.loc[missing,'smoker']=np.nan
        workers.loc[1,'age']=12;defects.append('W0002 age below eligibility range')
        if len(follow):follow.iloc[0,follow.columns.get_loc('time_years')]=-1;defects.append('one negative follow-up time')
        measurements.loc[0,'concentration_mg_m3']=-.1;measurements.loc[1,'unit']='ug/m3'
        workers=pd.concat([workers,workers.iloc[[0]]],ignore_index=True)
        defects+=['one exact duplicated worker row','one negative measurement','one invalid measurement unit','approximately three percent missing smoking']
    for name,df in [('workers',workers),('followup',follow),('measurements',measurements)]:df.to_csv(out/(name+'.csv'),index=False,float_format='%.12g')
    (out/'truth.json').write_text(json.dumps(dict(n_workers=n,seed=seed,synthetic=True,lag_years=5,administrative_years=8,weibull_shape=1.4,baseline_scale=.012,beta=[float(np.log(z)) for z in [1.7,1.3,1.5]],covariates=['log2_1plus_cumulative','age_per_decade','smoker'],hazard_ratios=[1.7,1.3,1.5],defects=defects),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=500);p.add_argument('--seed',type=int,default=2026);p.add_argument('--out',default='raw');a=p.parse_args();simulate(a.n,a.seed,a.out)
