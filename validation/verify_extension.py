"""Independent output integrity checks for the longitudinal scenario."""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from longitudinal import validate

root = Path(sys.argv[1])
p = root / 'longitudinal'
q = json.loads((p / 'qc.json').read_text())
s = json.loads((p / 'summary.json').read_text())
for name, expected in q['sha256'].items():
    assert hashlib.sha256((p / 'raw' / name).read_bytes()).hexdigest() == expected
i = pd.read_csv(p / 'raw/intervals.csv')
w = pd.read_csv(p / 'raw/subjects.csv')
h = pd.read_csv(p / 'raw/history.csv')
checked = validate(i, w, h)
for key in ['subjects', 'intervals', 'disease', 'deaths', 'censored']:
    assert checked[key] == s[key]
assert s['subjects'] == s['disease'] + s['deaths'] + s['censored']
assert np.isclose((i.stop-i.start).sum(), s['person_years'])
for _, g in pd.read_csv(p / 'aj_curves.csv').groupby('group'):
    assert np.allclose(g.survival + g.disease + g.death, 1)
    assert (np.diff(g.disease) >= -1e-12).all()
    assert (np.diff(g.death) >= -1e-12).all()
    assert (g.naive_one_minus_km >= g.disease-1e-12).all()
c = pd.read_csv(p / 'cif_intervals.csv')
assert ((c.lower >= 0) & (c.upper <= 1) & (c.lower <= c.upper)).all()
for m in s['models']:
    assert m['events'] == int(((i.cause == m['cause']) & i.smoker.notna()).sum())
    assert m['status'] == ('fitted' if m['events'] >= 30 else 'skipped')
    if m['status'] == 'fitted':
        assert m['score_max'] < 1e-4
        assert len(pd.read_csv(p / f"cox_cause{m['cause']}.csv")) == 3
for name in ['exposure_history', 'cumulative_incidence', 'cif_by_industry', 'time_varying_hr']:
    assert (root / 'longitudinal_figures' / f'{name}.png').stat().st_size > 1000
assert (root / 'report/OccuSurv_report.docx').stat().st_size > 10000
print(f"PASS {root}: hashes, interval QC, counts, person-time, probability mass, monotonicity, confidence bounds, model gates, figures and report")
