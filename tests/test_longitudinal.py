import unittest,tempfile,sys
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import approx_fprime
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import CumIncidenceRight
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'bin'))
from longitudinal import aj,lagged,validate,generate,fit
from simulate import simulate
from qc import validate as baseline_qc
class LongitudinalTests(unittest.TestCase):
 def fixture(self,t):
  t=Path(t);simulate(40,2026,t/'baseline');baseline_qc(t/'baseline',t/'clean');generate(t/'clean/workers.csv',2027,t/'raw');return [pd.read_csv(t/'raw'/name) for name in ['intervals.csv','subjects.csv','history.csv']]
 def test_aj_hand_calculation(self):
  d=aj([1,2,3],[1,2,0]);self.assertAlmostEqual(d.disease.iloc[-1],1/3);self.assertAlmostEqual(d.death.iloc[-1],1/3);self.assertAlmostEqual(d.survival.iloc[-1],1/3)
 def test_aj_matches_library_with_ties(self):
  t=[1,2,2,3,4,5];c=[1,2,1,0,2,0];a=aj(t,c);b=CumIncidenceRight(t,c)
  for j,col in enumerate(['disease','death']):np.testing.assert_allclose(a.set_index('time').loc[b.times,col],b.cinc[j])
 def test_mass_and_monotonicity(self):
  d=aj([1,2,3,4],[2,1,2,0]);np.testing.assert_allclose(d.survival+d.disease+d.death,1);self.assertTrue((np.diff(d.disease)>=0).all());self.assertGreater(d.naive_one_minus_km.iloc[-1],d.disease.iloc[-1])
 def test_all_censored(self):
  d=aj([1,2],[0,0]);self.assertEqual(d.disease.sum()+d.death.sum(),0)
 def test_lag_no_future_leak(self):
  e=[.1]*8;self.assertAlmostEqual(lagged(.2,10,e,0),1.);self.assertAlmostEqual(lagged(.2,10,e,6),2.1);e[7]=999;self.assertAlmostEqual(lagged(.2,10,e,6),2.1)
 def test_valid_histories(self):
  with tempfile.TemporaryDirectory() as t:
   i,s,h=self.fixture(t);self.assertEqual(validate(i,s,h)['status'],'PASS')
 def test_overlap_fails(self):
  with tempfile.TemporaryDirectory() as t:
   i,s,h=self.fixture(t);ix=i.index[i.start>0][0];i.loc[ix,'start']-=.1
   with self.assertRaises(ValueError):validate(i,s,h)
 def test_post_event_fails(self):
  with tempfile.TemporaryDirectory() as t:
   i,s,h=self.fixture(t);ix=i.groupby('worker_id').filter(lambda g:len(g)>1).index[0];i.loc[ix,'cause']=1
   with self.assertRaises(ValueError):validate(i,s,h)
 def test_exposure_assignment_fails(self):
  with tempfile.TemporaryDirectory() as t:
   i,s,h=self.fixture(t);i.loc[0,'log2_1plus_lagged']+=1
   with self.assertRaises(ValueError):validate(i,s,h)
 def test_terminal_mismatch_fails(self):
  with tempfile.TemporaryDirectory() as t:
   i,s,h=self.fixture(t);s.loc[0,'time']+=.1
   with self.assertRaises(ValueError):validate(i,s,h)
 def test_low_cause_events_skip(self):
  with tempfile.TemporaryDirectory() as t:
   i,s,h=self.fixture(t);self.assertEqual(fit(i,2,Path(t))['status'],'skipped')
 def test_start_stop_likelihood(self):
  # Events on and between interval boundaries; exactly one active row per worker.
  start=np.array([0.,1.,0.,2.,0.]);stop=np.array([1.,3.,2.,4.,2.]);e=np.array([0,1,0,1,1]);x=np.array([[0.,1.],[1.,1.],[2.,0.],[3.,0.],[.5,1.]])
  model=PHReg(stop,x,status=e,entry=np.nextafter(start,np.inf));b=np.array([.2,-.3])
  def direct(beta):return sum(x[i]@beta-np.log(np.exp(x[(start<stop[i])&(stop>=stop[i])]@beta).sum()) for i in range(len(stop)) if e[i])
  self.assertAlmostEqual(model.loglike(b),direct(b),places=9);np.testing.assert_allclose(model.score(b),approx_fprime(b,direct,1e-7),atol=1e-6)
if __name__=='__main__':unittest.main()
