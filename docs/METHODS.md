# Statistical methodology

All parameters below describe a synthetic data-generating mechanism. No NIOH observations are used.

## Exposure

$$
C_i=c_i\max(d_i-5,0),\quad x_i=\log_2(1+C_i/C_0),\quad C_0=1\ {\rm mg\,year/m^3}
$$

Historical concentration c is constant over employment duration d in this simulation. Five years of recent exposure are excluded from the baseline cumulative metric. The reference unit C0 makes the logarithm dimensionless.

## Hazard

$$
h_i(t)=1.4(0.012)t^{0.4}\exp(\eta_i),\quad \eta_i=\log(1.7)x_i+\log(1.3)a_i+\log(1.5)s_i
$$

Time is measured in years; a=(age-45)/10 and s is smoking status. All constants are illustrative generating assumptions, not estimates for South African workers.

## Event

$$
T_i=\left[\frac{-\log U_i}{0.012\exp(\eta_i)}\right]^{1/1.4},\quad Y_i=\min(T_i,C_i^*),\quad \delta_i=I(T_i\leq C_i^*)
$$

U is uniform on (0,1); independent censoring C* is the minimum of eight years and an exponential time with mean 25 years. C* is distinct from cumulative exposure C. Prevalent cases do not enter the incident cohort.

## Wilson

$$
\frac{\hat p+z^2/(2n)\ \pm\ z\sqrt{\hat p(1-\hat p)/n+z^2/(4n^2)}}{1+z^2/n},\quad \hat p=k/n,\quad z=1.96
$$

The baseline prevalence proportion uses a Wilson 95 percent interval. It describes the eligible cross-sectional sample and is not an incidence estimate.

## Rate

$$
\hat r=1000D/PY,\quad CI=\frac{1000}{2PY}\left[\chi^2_{2D,0.025},\chi^2_{2(D+1),0.975}\right]
$$

D is the number of incident events and PY is observed person-time. The lower endpoint is zero for D=0. These are conventional Poisson count intervals conditional on person-time; heterogeneity and workplace dependence can limit their coverage.

## Km

$$
\hat S(t)=\prod_{t_j\leq t}(1-d_j/n_j),\quad G(t)=\sum_{t_j\leq t}\frac{d_j}{n_j(n_j-d_j)}
$$

At each distinct event time, nj is the risk set and dj the number of diagnoses. Greenwood G supports log-minus-log pointwise confidence intervals: transform log(-log S), use standard error sqrt(G)/abs(log S), then back-transform. Boundary survival values use boundary limits.

## Cox

$$
\ell(\beta)=\sum_{i:\delta_i=1}\left[x_i^T\beta-\log\left(\sum_{j:Y_j\geq Y_i}\exp(x_j^T\beta)\right)\right]
$$

The Cox partial log likelihood uses the exposure transform, age per decade and smoking. Continuous generated event times have no designed ties; Breslow handling is configured for ties. No intercept is fitted because the baseline hazard absorbs it.

## Newton

$$
\beta_{m+1}=\beta_m+I(\beta_m)^{-1}U(\beta_m),\quad V=I(\hat\beta)^{-1}
$$

Statsmodels maximises the partial likelihood by Newton optimisation, using score U and observed information I. We require a finite solution and maximum absolute score below 1e-4. V is model-based covariance, conditional on the specified covariates and independent workers.

## Hr

$$
HR_k=\exp(\hat\beta_k),\quad CI_k=\exp\left(\hat\beta_k\pm1.96\sqrt{V_{kk}}\right)
$$

The exposure HR is per one-unit increase in log2(1+C/C0), not per doubling of raw exposure. Wald P values are exploratory and unadjusted. The model is skipped below 30 complete-case events, a computational safeguard rather than a sample-size justification.

## Residual

$$
r_i=x_i-\frac{\sum_{j\in R(Y_i)}x_j\exp(x_j^T\hat\beta)}{\sum_{j\in R(Y_i)}\exp(x_j^T\hat\beta)}
$$

Schoenfeld residuals at event times are plotted against log time with LOWESS smoothing. Spearman correlations are descriptive only. These diagnostics are not a formal global proportional-hazards test and do not establish that proportional hazards holds.

## Additional generating assumptions

Twelve workplaces are distributed across three industries. Workplace effects alter exposure only, not residual disease hazard. Concentration is lognormal with industry medians 0.04, 0.07 and 0.12 mg/m3, workplace log SD 0.18, and worker log SD 0.35. Age is uniform from 25 to 65; historical employment is the smaller of Uniform(3,30) years and age minus 18. Three baseline measurements per worker have multiplicative lognormal error with log SD 0.25. They illustrate monitoring and do not replace the known simulated historical concentration in the outcome model.

Smoking is Bernoulli(0.35). Baseline disease is Bernoulli with probability expit(-3.4 + 0.5 x + 0.3 a). Follow-up disease times are drawn independently of baseline disease conditional on covariates. This specifies a conditional cohort model; it does not reconstruct lifetime silicosis onset. Three percent missing smoking is introduced completely at random, plus deliberate QC defects. Missing smoking is retained in descriptive analyses and excluded from regression; no imputation is performed.

There is no competing mortality, migration-related selection, unmeasured shared workplace frailty or treatment intervention. These restrictions are essential to interpreting the demonstration. A single small simulation checks computation, not estimator bias or interval coverage.

## Score and observed information

For a risk set at event time t, define weighted sums for k=0,1,2:

$$
S^{(0)}(t,\beta)=\sum_{j\in R(t)}e^{X_j^T\beta},\quad
S^{(1)}(t,\beta)=\sum_{j\in R(t)}X_j e^{X_j^T\beta},\quad
S^{(2)}(t,\beta)=\sum_{j\in R(t)}X_jX_j^T e^{X_j^T\beta}.
$$

For the continuous event times used here, the score and observed information are:

$$
U(\beta)=\sum_{i:\delta_i=1}\left[X_i-\frac{S^{(1)}(Y_i,\beta)}{S^{(0)}(Y_i,\beta)}\right],
$$

$$
I(\beta)=\sum_{i:\delta_i=1}\left[\frac{S^{(2)}(Y_i,\beta)}{S^{(0)}(Y_i,\beta)}-\bar X(Y_i,\beta)\bar X(Y_i,\beta)^T\right],
\qquad \bar X=S^{(1)}/S^{(0)}.
$$

Thus U is the gradient of the partial log likelihood and I is its negative Hessian. The unit tests compare the library log likelihood against an explicit risk-set sum and its score against a finite-difference gradient. They also compare the implemented Kaplan–Meier estimates against Statsmodels on a small dataset with tied times. These numerical checks do not establish statistical coverage.
