#!/usr/bin/env python3
"""Separate longitudinal scenario on the baseline eligible workers."""
import argparse,json,hashlib,warnings
from pathlib import Path
import numpy as np,pandas as pd
from statsmodels.duration.hazard_regression import PHReg

TERMS=['log2_1plus_lagged','age_per_decade','smoker']
TRUE={1:np.log([1.7,1.3,1.5]),2:np.log([1.,1.5,1.2])}

def lagged(c,d,exposure,k,lag=5):
    # Integral of exposure from employment entry to k-lag; held at annual boundaries.
    historical=c*max(d-max(lag-k,0),0)
    return historical+sum(exposure[:max(k-lag,0)])

def generate(workers,seed,out):
    out=Path(out);out.mkdir(exist_ok=True);rng=np.random.default_rng(seed)
    w=pd.read_csv(workers);w=w.loc[w.baseline_disease==0].sort_values('worker_id');intervals=[];history=[];subjects=[]
    for row in w.itertuples():
        latent_smoke=int(row.smoker) if pd.notna(row.smoker) else int(rng.binomial(1,.35))
        exposure=row.historical_mg_m3*np.exp(rng.normal(0,.25,8))*np.array([1,1,1,.65,.65,.65,.65,.65])
        censor=min(8.,rng.exponential(25));observed=censor;cause=0
        for k in range(8):
            history.append(dict(worker_id=row.worker_id,industry=row.industry,year=k,concentration=exposure[k],lagged_cumulative=lagged(row.historical_mg_m3,row.employment_years,exposure,k)))
            if k>=observed:continue
            x=np.array([np.log2(1+lagged(row.historical_mg_m3,row.employment_years,exposure,k)),row.age_per_decade,latent_smoke])
            rates=np.array([.018*np.exp(TRUE[1]@x),.025*np.exp(TRUE[2]@x)])
            wait=rng.exponential(1/rates.sum());end=min(k+1.,censor)
            if k+wait<end:
                observed=k+wait;cause=1 if rng.uniform()<rates[0]/rates.sum() else 2;end=observed
            intervals.append(dict(worker_id=row.worker_id,industry=row.industry,start=float(k),stop=end,cause=cause if end==observed else 0,log2_1plus_lagged=x[0],age_per_decade=x[1],smoker=row.smoker))
        subjects.append(dict(worker_id=row.worker_id,industry=row.industry,time=observed,cause=cause,smoker=row.smoker))
    for name,data in [('intervals',intervals),('history',history),('subjects',subjects)]:pd.DataFrame(data).to_csv(out/(name+'.csv'),index=False)
    (out/'truth.json').write_text(json.dumps(dict(seed=seed,n=len(w),lag_years=5,years=8,disease_baseline=.018,death_baseline=.025,hr_disease=[1.7,1.3,1.5],hr_death=[1.,1.5,1.2],terms=TERMS,bootstrap_replicates=200,scenario='new outcomes on same baseline disease-free workers',planned_concentrations='latent complete schedules; post-event values never used in risk sets'),indent=2))

def validate(intervals,subjects,history):
    if subjects.worker_id.isna().any() or subjects.worker_id.duplicated().any():raise ValueError('Invalid subject IDs')
    if set(intervals.worker_id)!=set(subjects.worker_id) or set(history.worker_id)!=set(subjects.worker_id):raise ValueError('Unmatched subject histories')
    for table,cols in [(intervals,['start','stop','cause','log2_1plus_lagged','age_per_decade']),(subjects,['time','cause']),(history,['year','concentration','lagged_cumulative'])]:
        if not np.isfinite(table[cols].to_numpy()).all():raise ValueError('Nonfinite required values')
    if not intervals.cause.isin([0,1,2]).all() or not subjects.cause.isin([0,1,2]).all():raise ValueError('Invalid event code')
    if not (intervals.start.ge(0)&intervals.stop.gt(intervals.start)&intervals.stop.le(8)).all():raise ValueError('Invalid interval')
    if not (subjects.time.gt(0)&subjects.time.le(8)).all():raise ValueError('Invalid terminal time')
    if not (intervals.smoker.isna()|intervals.smoker.isin([0,1])).all():raise ValueError('Invalid smoking')
    if not (history.concentration.ge(0)&history.lagged_cumulative.ge(0)).all():raise ValueError('Negative exposure')
    ss=subjects.set_index('worker_id')
    for ident,g in intervals.groupby('worker_id'):
        g=g.sort_values('start');s=ss.loc[ident]
        if g.start.iloc[0]!=0 or not np.allclose(g.start.to_numpy()[1:],g.stop.to_numpy()[:-1],rtol=0,atol=1e-10):raise ValueError('Gap or overlap')
        if not np.isclose(g.stop.iloc[-1],s.time) or int(g.cause.iloc[-1])!=int(s.cause):raise ValueError('Terminal mismatch')
        if (g.cause.iloc[:-1]!=0).any():raise ValueError('Observation after event')
        if not (g.industry==s.industry).all() or g.age_per_decade.nunique()!=1:raise ValueError('Inconsistent subject covariates')
        h=history.loc[history.worker_id==ident].sort_values('year')
        if not np.array_equal(h.year.to_numpy(),np.arange(8)):raise ValueError('Incomplete or duplicated annual history')
        if np.any(np.diff(h.lagged_cumulative)<-1e-10):raise ValueError('Decreasing cumulative exposure')
        expected=np.log2(1+h.set_index('year').loc[g.start.astype(int),'lagged_cumulative'].to_numpy())
        if not np.allclose(g.log2_1plus_lagged,expected):raise ValueError('Future or inconsistent exposure assignment')
    return dict(status='PASS',subjects=len(subjects),intervals=len(intervals),disease=int((subjects.cause==1).sum()),deaths=int((subjects.cause==2).sum()),censored=int((subjects.cause==0).sum()),person_years=float(subjects.time.sum()),model_complete_subjects=int(subjects.smoker.notna().sum()),checks=['IDs','finite values','event codes','positive nonoverlapping contiguous intervals','terminal first-event matching','annual history','monotone exposure','covariate schedule agreement'])

def aj(time,cause):
    t=np.asarray(time);c=np.asarray(cause);times,ix,counts=np.unique(t,return_inverse=True,return_counts=True)
    risk=len(t)-np.cumsum(counts)+counts
    d1=np.bincount(ix,weights=c==1,minlength=len(times));d2=np.bincount(ix,weights=c==2,minlength=len(times))
    survival=np.cumprod(1-(d1+d2)/risk);before=np.r_[1.,survival[:-1]]
    f1=np.cumsum(before*d1/risk);f2=np.cumsum(before*d2/risk)
    net=1-np.cumprod(1-d1/risk)
    return pd.DataFrame(dict(time=np.r_[0,times],survival=np.r_[1,survival],disease=np.r_[0,f1],death=np.r_[0,f2],naive_one_minus_km=np.r_[0,net]))

def at_grid(curve,grid,column):
    ix=np.searchsorted(curve.time,grid,side='right')-1
    return curve[column].to_numpy()[np.maximum(ix,0)]

def fit(intervals,cause,out):
    d=intervals.dropna(subset=TERMS);e=(d.cause==cause).astype(int);events=int(e.sum());n=d.worker_id.nunique()
    result=dict(cause=cause,events=events,subjects=n,intervals=len(d),status='skipped',reason='fewer than 30 cause-specific events')
    if events<30:return result
    # At a boundary start==t, exclude the entering interval under (start,stop].
    entry=np.nextafter(d.start.to_numpy(),np.inf)
    model=PHReg(d.stop.to_numpy(),d[TERMS].to_numpy(),status=e.to_numpy(),entry=entry,ties='breslow')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always');res=model.fit(maxiter=100,disp=False)
    score=float(np.max(np.abs(model.score(res.params))))
    if not np.isfinite(res.params).all() or not np.isfinite(res.bse).all() or score>1e-4:raise ValueError('Cause-specific fit failed convergence')
    ci=res.conf_int();rows=[dict(cause=cause,term=term,beta=res.params[k],se=res.bse[k],hr=np.exp(res.params[k]),lower=np.exp(ci[k,0]),upper=np.exp(ci[k,1]),p=res.pvalues[k],assumed_hr=np.exp(TRUE[cause][k])) for k,term in enumerate(TERMS)]
    pd.DataFrame(rows).to_csv(out/f'cox_cause{cause}.csv',index=False)
    return dict(cause=cause,events=events,subjects=n,intervals=len(d),status='fitted',score_max=score,warnings=[str(x.message) for x in caught],covariance='model-based counting-process information; independent workers',ties='breslow')

def run(workers,seed,out):
    out=Path(out);out.mkdir(exist_ok=True);raw=out/'raw';generate(workers,seed,raw)
    i=pd.read_csv(raw/'intervals.csv');s=pd.read_csv(raw/'subjects.csv');h=pd.read_csv(raw/'history.csv');qc=validate(i,s,h)
    qc['sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in raw.iterdir()};(out/'qc.json').write_text(json.dumps(qc,indent=2))
    curves=[];endpoints=[];rng=np.random.default_rng(seed+1);grid=np.linspace(0,8,81);B=200
    for name,g in [('Overall',s)]+list(s.groupby('industry')):
        curve=aj(g.time,g.cause);curve['group']=name;curves.append(curve)
        draws={j:[] for j in ['disease','death']}
        for b in range(B):
            sampled=g.iloc[rng.integers(0,len(g),len(g))];bc=aj(sampled.time,sampled.cause)
            for j in draws:draws[j].append(at_grid(bc,grid,j))
        for j in draws:
            low,high=np.quantile(np.array(draws[j]),[.025,.975],axis=0);values=at_grid(curve,grid,j)
            for t,v,l,u in zip(grid,values,low,high):endpoints.append(dict(group=name,cause=j,time=t,cif=v,lower=l,upper=u,bootstrap_replicates=B))
    pd.concat(curves).to_csv(out/'aj_curves.csv',index=False);pd.DataFrame(endpoints).to_csv(out/'cif_intervals.csv',index=False)
    statuses=[fit(i,j,out) for j in [1,2]]
    result=dict(scenario='longitudinal competing events',seed=seed,**{k:qc[k] for k in ['subjects','intervals','disease','deaths','censored','person_years','model_complete_subjects']},models=statuses,bootstrap_replicates=B,worker_sampling='iid within each displayed group; no residual workplace frailty assumed')
    (out/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('workers');p.add_argument('--seed',type=int,default=2027);p.add_argument('--out',default='longitudinal');a=p.parse_args();run(a.workers,a.seed,a.out)
