import argparse,json
from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('--out',default='longitudinal_figures');a=p.parse_args();src=Path(a.source);out=Path(a.out);out.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
h=pd.read_csv(src/'raw/history.csv');sub=pd.read_csv(src/'raw/subjects.csv');curves=pd.read_csv(src/'aj_curves.csv');ci=pd.read_csv(src/'cif_intervals.csv');colors=['#176B87','#D07C32','#41876B']
# Only show measurements scheduled while workers remain under observation.
h=h.merge(sub[['worker_id','time']],on='worker_id',validate='many_to_one');h=h.loc[h.year<h.time]
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for (name,g),color in zip(h.groupby('industry'),colors):
 for ax,col in zip(axes,['concentration','lagged_cumulative']):
  med=g.groupby('year')[col].median();ax.step(med.index,med,where='post',c=color,label=name)
axes[0].set(xlabel='Year since baseline',ylabel='Concentration (mg/m3)',title='Annual exposure among observed workers');axes[1].set(xlabel='Year since baseline',ylabel='Lagged cumulative exposure (mg year/m3)',title='Five-year lag at annual boundaries');axes[0].legend(fontsize=8);fig.savefig(out/'exposure_history.png');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4.2),layout='constrained')
for cause,color in zip(['disease','death'],colors):
 d=ci.loc[(ci.group=='Overall')&(ci.cause==cause)];ax.plot(d.time,d.cif,label=cause.capitalize(),c=color);ax.fill_between(d.time,d.lower,d.upper,color=color,alpha=.15)
c=curves.loc[curves.group=='Overall'];ax.step(c.time,c.naive_one_minus_km,where='post',ls='--',c='#666666',label='Disease: 1 - KM with deaths censored')
ax.set(xlabel='Years since baseline',ylabel='Cumulative probability',ylim=(0,.65),title='Aalen Johansen cumulative incidence');ax.legend(fontsize=8);fig.savefig(out/'cumulative_incidence.png');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
for (name,g),color in zip(ci.loc[(ci.group!='Overall')&(ci.cause=='disease')].groupby('group'),colors):
 ax.plot(g.time,g.cif,c=color,label=name);ax.fill_between(g.time,g.lower,g.upper,color=color,alpha=.1)
ax.set(xlabel='Years since baseline',ylabel='Disease cumulative incidence',title='Competing-risk disease probability by industry',ylim=(0,.65));ax.legend(fontsize=8);fig.savefig(out/'cif_by_industry.png');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained');summary=json.loads((src/'summary.json').read_text())
for ax,j,label in zip(axes,[1,2],['Diagnosis','Competing death']):
 f=src/f'cox_cause{j}.csv'
 if f.exists():
  d=pd.read_csv(f);y=np.arange(3);ax.errorbar(d.hr,y,xerr=[d.hr-d.lower,d.upper-d.hr],fmt='o',capsize=4,c=colors[0]);ax.scatter(d.assumed_hr,y,c=colors[1],marker='x',label='Assumed HR');ax.axvline(1,c='grey',ls='--');ax.set(xscale='log',yticks=y,yticklabels=['Lagged exposure transform','Age per decade','Smoker'],xlabel='Cause-specific HR',title=label);ax.legend(fontsize=8)
 else:ax.text(.5,.5,'Model skipped\nFewer than 30 events',ha='center',transform=ax.transAxes);ax.set_title(label);ax.set_axis_off()
fig.savefig(out/'time_varying_hr.png');plt.close(fig)
