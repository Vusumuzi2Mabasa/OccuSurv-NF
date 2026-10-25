#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
p=argparse.ArgumentParser();p.add_argument('analysis');p.add_argument('--out',default='figures');a=p.parse_args();src=Path(a.analysis);out=Path(a.out);out.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
colors=['#176B87','#D07C32','#41876B']
m=pd.read_csv(src/'exposure_means.csv');names=sorted(m.industry.unique())
fig,ax=plt.subplots(figsize=(8,3.7),layout='constrained');data=[m.loc[m.industry==s,'mean_measured_mg_m3'] for s in names]
ax.boxplot(data,tick_labels=names,showfliers=False);rng=np.random.default_rng(19)
for i,(d,color) in enumerate(zip(data,colors)):ax.scatter(i+1+rng.uniform(-.14,.14,len(d)),d,c=color,s=8,alpha=.35)
ax.set(ylabel='Worker mean silica concentration (mg/m3)',title='Synthetic exposure measurements');fig.savefig(out/'exposure.png');plt.close(fig)
r=pd.read_csv(src/'incidence.csv').iloc[1:];fig,ax=plt.subplots(figsize=(8,3.5),layout='constrained');ax.errorbar(r.rate_per_1000,np.arange(len(r)),xerr=[r.rate_per_1000-r.lower,r.upper-r.rate_per_1000],fmt='o',capsize=4,c=colors[0]);ax.set(yticks=np.arange(len(r)),yticklabels=r.industry,xlabel='Incident diagnoses per 1000 person-years',title='Crude rates with 95% Poisson intervals');fig.savefig(out/'incidence.png');plt.close(fig)
k=pd.read_csv(src/'km.csv');fig,ax=plt.subplots(figsize=(8,3.7),layout='constrained')
for name,color in zip(names,colors):
 d=k.loc[k.industry==name];ax.step(d.time,d.survival,where='post',label=name,c=color);ax.fill_between(d.time,d.lower,d.upper,step='post',color=color,alpha=.1)
ax.set(xlabel='Years since baseline',ylabel='Estimated disease-free survival',ylim=(0,1.02),xlim=(0,8),title='Kaplan Meier estimates with pointwise 95% intervals');ax.legend(fontsize=8,loc='lower left');fig.savefig(out/'survival.png');plt.close(fig)
s=json.loads((src/'analysis.json').read_text())
if s['model']['status']=='fitted':
 d=pd.read_csv(src/'cox.csv');fig,ax=plt.subplots(figsize=(8,3.5),layout='constrained');y=np.arange(len(d));ax.errorbar(d.hr,y,xerr=[d.hr-d.lower,d.upper-d.hr],fmt='o',label='Estimate and 95% Wald interval',c=colors[0],capsize=4);ax.scatter(d.assumed_hr,y,marker='x',s=65,c=colors[1],label='Assumed generating value');ax.axvline(1,c='grey',ls='--');ax.set(xscale='log',yticks=y,yticklabels=['Historical exposure transform','Age per 10 years','Smoker versus nonsmoker'],xlabel='Adjusted hazard ratio',title='Cox proportional hazards model');ax.legend(fontsize=8);fig.savefig(out/'cox.png');plt.close(fig)
 r=pd.read_csv(src/'schoenfeld.csv');fig,axes=plt.subplots(1,3,figsize=(10,3.3),layout='constrained')
 for ax,(term,d) in zip(axes,r.groupby('term')):
  x=np.log(d.time.to_numpy());y=d.residual.to_numpy();ax.scatter(x,y,s=8,alpha=.45,c=colors[0]);fit=lowess(y,x,frac=.65);ax.plot(fit[:,0],fit[:,1],c=colors[1]);ax.axhline(0,c='grey',lw=.7);ax.set(title=term,xlabel='Log event time',ylabel='Schoenfeld residual')
 fig.savefig(out/'ph_diagnostics.png');plt.close(fig)
(out/'versions.json').write_text(json.dumps(dict(matplotlib=matplotlib.__version__),indent=2))
