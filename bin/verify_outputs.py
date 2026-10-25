#!/usr/bin/env python3
import argparse,csv,hashlib,json,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('results',type=Path);a=p.parse_args();r=a.results
s=json.loads((r/'analysis/analysis.json').read_text());q=json.loads((r/'clean/qc.json').read_text())
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
for name,digest in q['input_sha256'].items():assert hashlib.sha256((r/'raw'/name).read_bytes()).hexdigest()==digest,name
w=read(r/'clean/workers.csv');c=read(r/'clean/cohort.csv');rates=read(r/'analysis/incidence.csv')
assert len(w)==s['workers'] and len(c)==s['incident_cohort']
assert all(float(x['baseline_disease'])==0 and 0<float(x['time_years'])<=8 for x in c)
assert sum(int(x['event']) for x in c)==s['events']
assert math.isclose(sum(float(x['time_years']) for x in c),s['person_years'],rel_tol=1e-10)
assert math.isclose(float(rates[0]['rate_per_1000']),1000*s['events']/s['person_years'],rel_tol=1e-10)
if s['model']['status']=='fitted':
 assert len(read(r/'analysis/cox.csv'))==3 and s['model']['events']>=30 and s['model']['score_max']<1e-4
 expected=['exposure','incidence','km','cox_forest','ph_diagnostics']
else:
 assert not (r/'analysis/cox.csv').exists();expected=['exposure','incidence','km']
for name in expected:assert (r/'figures'/(name+'.png')).stat().st_size>1000
assert (r/'report/OccuSurv_report.docx').stat().st_size>1000
print(json.dumps(dict(status='PASS',workers=s['workers'],events=s['events'],model=s['model']['status'],checks=['raw hashes unchanged','cohort eligibility','event counts','person time','rates','model gate','figures','Word report']),indent=2))
