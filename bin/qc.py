#!/usr/bin/env python3
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd

def validate(source,out):
    source=Path(source);out=Path(out);out.mkdir(exist_ok=True);audit=[]
    def load(name,cols):
        d=pd.read_csv(source/(name+'.csv'))
        if not set(cols)<=set(d.columns):raise ValueError('Missing columns in '+name)
        return d
    def issue(table,idx,identifier,reason,action='exclude'):
        audit.append(dict(table=table,source_row=int(idx)+2,identifier=str(identifier),reason=reason,action=action))
    w=load('workers',['worker_id','workplace','industry','age','smoker','employment_years','historical_mg_m3','baseline_disease'])
    duplicate=w.duplicated();
    for i in w.index[duplicate]:issue('workers',i,w.loc[i,'worker_id'],'exact_duplicate')
    w=w.loc[~duplicate].copy()
    if w.worker_id.isna().any() or w.worker_id.duplicated().any():raise ValueError('Missing or conflicting worker IDs')
    valid=pd.Series(True,index=w.index)
    checks={'invalid_age':w.age.between(18,80),'invalid_employment':(w.employment_years>=0)&(w.employment_years<=w.age-18+1e-8),'invalid_history':np.isfinite(w.historical_mg_m3)&(w.historical_mg_m3>=0),'invalid_baseline':w.baseline_disease.isin([0,1]),'invalid_smoking':w.smoker.isna()|w.smoker.isin([0,1]),'invalid_group':w.workplace.notna()&w.industry.isin(['Construction','Brick manufacture','Stone fabrication'])}
    for reason,ok in checks.items():
        for i in w.index[~ok]:issue('workers',i,w.loc[i,'worker_id'],reason)
        valid &=ok
    for i in w.index[w.smoker.isna()]:issue('workers',i,w.loc[i,'worker_id'],'missing_smoking','retain_descriptive_exclude_model')
    w=w.loc[valid].copy();w['cumulative_mg_year_m3']=w.historical_mg_m3*np.maximum(w.employment_years-5,0);w['log2_1plus_cumulative']=np.log2(1+w.cumulative_mg_year_m3);w['age_per_decade']=(w.age-45)/10
    f=load('followup',['worker_id','time_years','event'])
    if f.worker_id.duplicated().any():raise ValueError('Duplicated follow-up IDs')
    eligible=set(w.loc[w.baseline_disease==0,'worker_id']);good=pd.Series(True,index=f.index)
    for reason,ok in {'not_eligible':f.worker_id.isin(eligible),'invalid_time':np.isfinite(f.time_years)&f.time_years.gt(0)&f.time_years.le(8),'invalid_event':f.event.isin([0,1])}.items():
        for i in f.index[~ok]:issue('followup',i,f.loc[i,'worker_id'],reason)
        good &=ok
    f=f.loc[good].copy()
    for identifier in sorted(eligible-set(f.worker_id)):issue('followup',-1,identifier,'no_valid_followup','exclude_incidence')
    c=w.merge(f,on='worker_id',how='inner',validate='one_to_one');c['model_complete']=c.smoker.notna()
    m=load('measurements',['worker_id','sample_id','concentration_mg_m3','duration_hours','unit'])
    if m.sample_id.isna().any() or m.sample_id.duplicated().any():raise ValueError('Invalid measurement IDs')
    good=pd.Series(True,index=m.index)
    for reason,ok in {'unknown_worker':m.worker_id.isin(w.worker_id),'invalid_concentration':np.isfinite(m.concentration_mg_m3)&m.concentration_mg_m3.ge(0),'invalid_duration':np.isfinite(m.duration_hours)&m.duration_hours.gt(0)&m.duration_hours.le(24),'invalid_unit':m.unit.eq('mg/m3')}.items():
        for i in m.index[~ok]:issue('measurements',i,m.loc[i,'sample_id'],reason)
        good &=ok
    m=m.loc[good].copy()
    for name,df in [('workers',w),('cohort',c),('measurements',m)]:df.to_csv(out/(name+'.csv'),index=False)
    pd.DataFrame(audit,columns=['table','source_row','identifier','reason','action']).to_csv(out/'audit.csv',index=False)
    stats=dict(workers_raw=len(pd.read_csv(source/'workers.csv')),workers_retained=len(w),followup_retained=len(c),model_complete=int(c.model_complete.sum()),measurements_raw=len(pd.read_csv(source/'measurements.csv')),measurements_retained=len(m),audit_entries=len(audit),input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.glob('*') if p.is_file()})
    (out/'qc.json').write_text(json.dumps(stats,indent=2));return stats
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('--out',default='clean');a=p.parse_args();validate(a.source,a.out)
