"""Fail the Nextflow task when R and independent Python estimators disagree."""
import sys,json,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from longitudinal import aj,at_grid,fit
from analyse import analyse

clean,raw,analysis,longitudinal=sys.argv[1:]
out=Path('crosscheck');out.mkdir(exist_ok=True)
ref=out/'python';ref.mkdir(exist_ok=True)
analyse(clean,Path(raw)/'truth.json',ref/'analysis')
lp=Path(longitudinal);s=pd.read_csv(lp/'raw/subjects.csv');i=pd.read_csv(lp/'raw/intervals.csv')
(ref/'longitudinal').mkdir(exist_ok=True)
statuses=[fit(i,k,ref/'longitudinal') for k in (1,2)]
rows=[]
def compare(label,a,b,tol=1e-7):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    assert a.shape==b.shape,(label,a.shape,b.shape)
    assert np.array_equal(np.isnan(a),np.isnan(b)),label
    delta=float(np.nanmax(np.abs(a-b))) if np.isfinite(a).any() else 0.
    passed=np.allclose(a,b,atol=tol,rtol=0,equal_nan=True)
    rows.append(dict(check=label,max_abs_difference=delta,tolerance=tol,passed=bool(passed)))
def table(label,a,b,keys,cols=None,tol=1e-7):
    a=pd.read_csv(a).sort_values(keys).reset_index(drop=True);b=pd.read_csv(b).sort_values(keys).reset_index(drop=True)
    assert a[keys].equals(b[keys]),label
    if cols is None:cols=list(a.select_dtypes('number').columns.difference(keys))
    compare(label,a[cols],b[cols],tol)
for n in ['prevalence','incidence','table1']:
    table(n,Path(analysis)/f'{n}.csv',ref/'analysis'/f'{n}.csv',['industry'])
table('exposure means',Path(analysis)/'exposure_means.csv',ref/'analysis/exposure_means.csv',['worker_id','industry'])
table('KM estimates and risk sets',Path(analysis)/'km.csv',ref/'analysis/km.csv',['industry','time'],['survival','at_risk','events'])
# Recompute the historical bounds with the exact critical value used by R.
# Keep the original Python table intact for audit, rather than widening tolerance.
pk=pd.read_csv(ref/'analysis/km.csv');parts=[]
for _,g in pk.groupby('industry',sort=False):
    g=g.copy();risk=g.at_risk.to_numpy();events=g.events.to_numpy()
    green=np.cumsum(np.divide(events,risk*(risk-events),out=np.zeros(len(g)),where=risk>events))
    surv=g.survival.to_numpy();valid=(surv>0)&(surv<1)
    low=surv.copy();high=surv.copy();se=np.sqrt(green[valid])/np.abs(np.log(surv[valid]));clog=np.log(-np.log(surv[valid]))
    low[valid]=np.exp(-np.exp(clog+norm.ppf(.975)*se));high[valid]=np.exp(-np.exp(clog-norm.ppf(.975)*se))
    g['lower']=low;g['upper']=high;parts.append(g)
pd.concat(parts).to_csv(ref/'analysis/km_exact_normal.csv',index=False)
table('KM bounds with exact normal critical value',Path(analysis)/'km.csv',ref/'analysis/km_exact_normal.csv',['industry','time'],['lower','upper'],tol=1e-12)
table('risk table',Path(analysis)/'risk_table.csv',ref/'analysis/risk_table.csv',['industry','time'])
rs=json.loads((Path(analysis)/'analysis.json').read_text());ps=json.loads((ref/'analysis/analysis.json').read_text())
assert rs['model']['status']==ps['model']['status']
if rs['model']['status']=='fitted':
    table('baseline Cox',Path(analysis)/'cox.csv',ref/'analysis/cox.csv',['term'])
    table('Schoenfeld rank correlations',Path(analysis)/'ph_diagnostics.csv',ref/'analysis/ph_diagnostics.csv',['term'])
rstatuses=json.loads((lp/'summary.json').read_text())['models']
for k,st in enumerate(statuses,1):
    assert st['status']==rstatuses[k-1]['status']
    if st['status']=='fitted':table(f'cause {k} Cox',lp/f'cox_cause{k}.csv',ref/'longitudinal'/f'cox_cause{k}.csv',['term'])
ac=pd.read_csv(lp/'aj_curves.csv');ci=pd.read_csv(lp/'cif_intervals.csv');indices=json.loads((lp/'bootstrap_indices.json').read_text())
for group,d in [('Overall',s)]+list(s.groupby('industry')):
    curve=aj(d.time,d.cause);r=ac[ac.group==group].sort_values('time')
    compare(f'{group} Aalen Johansen',r[['time','survival','disease','death','naive_one_minus_km']],curve,1e-12)
    assert np.allclose(r.survival+r.disease+r.death,1)
    ix=np.asarray(indices[group],dtype=int)-1
    draws={'disease':[],'death':[]};grid=np.linspace(0,8,81)
    for b in range(ix.shape[1]):
        sampled=d.iloc[ix[:,b]];bc=aj(sampled.time,sampled.cause)
        for cause in draws:draws[cause].append(at_grid(bc,grid,cause))
    for cause in draws:
        rr=ci[(ci.group==group)&(ci.cause==cause)].sort_values('time')
        bounds=np.quantile(draws[cause],[.025,.975],axis=0).T
        compare(f'{group} {cause} bootstrap',rr[['lower','upper']],bounds,1e-12)
pd.DataFrame(rows).to_csv(out/'comparison.csv',index=False)
result=dict(status='PASS' if all(r['passed'] for r in rows) else 'FAIL',checks=len(rows),baseline_model=rs['model']['status'],cause_models=[x['status'] for x in rstatuses],bootstrap='identical R-generated worker indices recomputed by Python')
(out/'summary.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result));assert result['status']=='PASS',rows
