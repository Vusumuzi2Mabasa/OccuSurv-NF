import sys,tempfile,unittest,json
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from simulate import simulate
from qc import validate
from analyse import wilson,rate_ci,km,cox_fit
class Tests(unittest.TestCase):
 def test_seed_reproducibility(self):
  with tempfile.TemporaryDirectory() as t:
   a=Path(t)/'a';b=Path(t)/'b';simulate(40,2026,a);simulate(40,2026,b)
   for f in a.iterdir():self.assertEqual(f.read_bytes(),(b/f.name).read_bytes())
 def test_zero_events_interval(self):
  r,lo,hi=rate_ci(0,100);self.assertEqual(r,0);self.assertEqual(lo,0);self.assertAlmostEqual(hi,36.8887945,places=5)
 def test_wilson_boundaries(self):
  for k in [0,10]:
   p,l,u=wilson(k,10);self.assertGreaterEqual(l,-1e-12);self.assertLessEqual(u,1+1e-12);self.assertLessEqual(l,p+1e-12);self.assertGreaterEqual(u,p-1e-12)
 def test_km_known(self):
  r=km([1,2,3],[1,0,1]);np.testing.assert_allclose(r.survival,[1,2/3,2/3,0])
 def test_km_all_censored(self):self.assertTrue((km([1,2],[0,0]).survival==1).all())
 def test_nonpositive_time_denominator(self):
  with self.assertRaises(ValueError):rate_ci(1,0)
 def test_qc_defects_and_hashes(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t)/'raw';c=Path(t)/'clean';simulate(40,2026,r);s=validate(r,c);self.assertEqual(s['workers_retained'],39)
   a=pd.read_csv(c/'audit.csv');self.assertTrue({'exact_duplicate','invalid_age','invalid_unit','invalid_concentration','invalid_time'}<=set(a.reason));self.assertEqual(len(s['input_sha256']),4)
 def test_conflicting_ids_fail(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t)/'raw';simulate(40,1,r,False);d=pd.read_csv(r/'workers.csv');d.loc[1,'worker_id']=d.loc[0,'worker_id'];d.to_csv(r/'workers.csv',index=False)
   with self.assertRaises(ValueError):validate(r,Path(t)/'clean')
 def test_missing_column_fail(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t)/'raw';simulate(40,1,r,False);d=pd.read_csv(r/'workers.csv');d.drop(columns='age').to_csv(r/'workers.csv',index=False)
   with self.assertRaises(ValueError):validate(r,Path(t)/'clean')
 def test_low_event_model_skips(self):
  c=pd.DataFrame(dict(log2_1plus_cumulative=[0,1],age_per_decade=[0,1],smoker=[0,1],event=[0,1],time_years=[2,1]))
  with tempfile.TemporaryDirectory() as t:self.assertEqual(cox_fit(c,Path(t))['status'],'skipped')
 def test_cumulative_exposure_formula(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t)/'raw';c=Path(t)/'clean';simulate(40,7,r,False);validate(r,c);d=pd.read_csv(c/'workers.csv');np.testing.assert_allclose(d.cumulative_mg_year_m3,d.historical_mg_m3*np.maximum(d.employment_years-5,0))
 def test_incident_cohort_excludes_prevalent(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t)/'raw';c=Path(t)/'clean';simulate(100,7,r);validate(r,c);d=pd.read_csv(c/'cohort.csv');self.assertTrue((d.baseline_disease==0).all());self.assertTrue((d.time_years>0).all())
 def test_cox_likelihood_and_score(self):
  from statsmodels.duration.hazard_regression import PHReg
  from scipy.optimize import approx_fprime
  t=np.array([1.,2.,3.,4.,5.]);e=np.array([1,0,1,1,0]);x=np.array([[0.,1.],[1.,0.],[2.,1.],[-1.,0.],[.5,1.]])
  model=PHReg(t,x,status=e);b=np.array([.2,-.3])
  def direct(beta):return sum(x[i]@beta-np.log(np.exp(x[t>=t[i]]@beta).sum()) for i in range(len(t)) if e[i])
  self.assertAlmostEqual(model.loglike(b),direct(b),places=10)
  np.testing.assert_allclose(model.score(b),approx_fprime(b,direct,1e-7),atol=1e-6)
 def test_km_matches_library(self):
  from statsmodels.duration.survfunc import SurvfuncRight
  t=[1,2,2,4,6];e=[1,0,1,1,0];a=km(t,e);b=SurvfuncRight(t,e)
  np.testing.assert_allclose(a.loc[a.events>0,'survival'],b.surv_prob)
if __name__=='__main__':unittest.main()
