#!/usr/bin/env python3
import argparse,json,platform,warnings
from pathlib import Path
import numpy as np,pandas as pd,scipy,statsmodels
from scipy.stats import chi2,norm,spearmanr
from statsmodels.duration.hazard_regression import PHReg

def wilson(k,n):
    z=norm.ppf(.975);p=k/n;den=1+z*z/n;center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return p,center-half,center+half

def rate_ci(k,py):
    if py<=0:raise ValueError('Nonpositive person-time')
    return 1000*k/py,0 if k==0 else 500*chi2.ppf(.025,2*k)/py,500*chi2.ppf(.975,2*(k+1))/py

def km(time,event):
    time=np.asarray(time);event=np.asarray(event);survival=1.;greenwood=0.;rows=[dict(time=0.,survival=1.,lower=1.,upper=1.,at_risk=len(time),events=0)]
    for t in np.unique(time):
        risk=int((time>=t).sum());d=int(((time==t)&(event==1)).sum())
        if d:
            survival*=1-d/risk
            if risk>d:greenwood+=d/(risk*(risk-d))
        if 0<survival<1:
            se=np.sqrt(greenwood)/abs(np.log(survival));a=np.log(-np.log(survival));lo=np.exp(-np.exp(a+1.96*se));hi=np.exp(-np.exp(a-1.96*se))
        else:lo=hi=survival
        rows.append(dict(time=float(t),survival=survival,lower=lo,upper=hi,at_risk=risk,events=d))
    return pd.DataFrame(rows)

def cox_fit(c,out):
    columns=['log2_1plus_cumulative','age_per_decade','smoker'];d=c.dropna(subset=columns).copy();events=int(d.event.sum())
    status=dict(n=len(d),events=events,parameters=3,minimum_events=30,status='skipped',reason='fewer than 30 events')
    if events<30:return status
    x=d[columns].to_numpy();t=d.time_years.to_numpy();e=d.event.to_numpy()
    if np.linalg.matrix_rank(x)<3:status['reason']='rank deficient';return status
    model=PHReg(t,x,status=e,ties='breslow')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always');fit=model.fit(maxiter=100,disp=False)
    score=float(np.max(np.abs(model.score(fit.params))))
    if not np.all(np.isfinite(fit.params)) or not np.all(np.isfinite(fit.bse)) or score>1e-4:raise RuntimeError('Cox fit did not meet convergence checks')
    ci=fit.conf_int();true=np.log([1.7,1.3,1.5]);rows=[]
    for k,col in enumerate(columns):rows.append(dict(term=col,beta=fit.params[k],se=fit.bse[k],hr=np.exp(fit.params[k]),lower=np.exp(ci[k,0]),upper=np.exp(ci[k,1]),p=fit.pvalues[k],assumed_hr=np.exp(true[k])))
    pd.DataFrame(rows).to_csv(out/'cox.csv',index=False)
    residuals=fit.schoenfeld_residuals;rr=[];diagnostics=[]
    for k,col in enumerate(columns):
        valid=(e==1)&np.isfinite(residuals[:,k]);rho=float(spearmanr(np.log(t[valid]),residuals[valid,k]).statistic)
        diagnostics.append(dict(term=col,spearman_rho=rho,n_events=int(valid.sum())))
        rr.extend(dict(term=col,time=float(a),residual=float(b)) for a,b in zip(t[valid],residuals[valid,k]))
    pd.DataFrame(rr).to_csv(out/'schoenfeld.csv',index=False);pd.DataFrame(diagnostics).to_csv(out/'ph_diagnostics.csv',index=False)
    return dict(status='fitted',n=len(d),events=events,parameters=3,minimum_events=30,score_max=score,ties='breslow',warnings=[str(x.message) for x in caught],covariance='model-based conditional independence',log_likelihood=float(fit.llf))

def analyse(source,truth,out):
    source=Path(source);out=Path(out);out.mkdir(exist_ok=True)
    w=pd.read_csv(source/'workers.csv');c=pd.read_csv(source/'cohort.csv');m=pd.read_csv(source/'measurements.csv')
    prevalence=[];incidence=[];descriptive=[]
    for name in ['Overall']+sorted(w.industry.unique()):
        ww=w if name=='Overall' else w.loc[w.industry==name];cc=c if name=='Overall' else c.loc[c.industry==name]
        pp=wilson(int(ww.baseline_disease.sum()),len(ww));rr=rate_ci(int(cc.event.sum()),cc.time_years.sum())
        prevalence.append(dict(industry=name,n=len(ww),cases=int(ww.baseline_disease.sum()),prevalence=pp[0],lower=pp[1],upper=pp[2]))
        incidence.append(dict(industry=name,n=len(cc),events=int(cc.event.sum()),person_years=cc.time_years.sum(),rate_per_1000=rr[0],lower=rr[1],upper=rr[2]))
        descriptive.append(dict(industry=name,n=len(ww),age_mean=ww.age.mean(),age_sd=ww.age.std(),smoking_observed=int(ww.smoker.notna().sum()),smoking_missing=int(ww.smoker.isna().sum()),smokers=int(ww.smoker.fillna(0).sum()),cumulative_median=ww.cumulative_mg_year_m3.median()))
    for name,rows in [('prevalence',prevalence),('incidence',incidence),('table1',descriptive)]:pd.DataFrame(rows).to_csv(out/(name+'.csv'),index=False)
    means=m.groupby('worker_id').concentration_mg_m3.mean().rename('mean_measured_mg_m3').reset_index().merge(w[['worker_id','industry']],validate='one_to_one');means.to_csv(out/'exposure_means.csv',index=False)
    curves=[];risk=[]
    for name,d in c.groupby('industry'):
        curve=km(d.time_years,d.event);curve['industry']=name;curves.append(curve)
        for t in [0,2,4,6,8]:risk.append(dict(industry=name,time=t,at_risk=int((d.time_years>=t).sum())))
    pd.concat(curves).to_csv(out/'km.csv',index=False);pd.DataFrame(risk).to_csv(out/'risk_table.csv',index=False)
    status=cox_fit(c,out)
    summary=dict(synthetic=True,workers=len(w),prevalent_cases=int(w.baseline_disease.sum()),incident_cohort=len(c),events=int(c.event.sum()),person_years=float(c.time_years.sum()),model=status,truth=json.loads(Path(truth).read_text()),versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__))
    (out/'analysis.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('truth');p.add_argument('--out',default='analysis');a=p.parse_args();analyse(a.source,a.truth,a.out)
