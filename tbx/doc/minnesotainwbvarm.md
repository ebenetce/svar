# svar.minnesotainwbvarm

Independent Normal-Wishart Minnesota prior for Bayesian VAR models.

The `svar.minnesotainwbvarm` model object specifies a Minnesota prior using the
`semiconjugatebvarm` parameterization. It supports a free cross-variable
tightness hyperparameter, $$\lambda_2$$, at the cost of losing the analytic
marginal likelihood available in the conjugate MNIW variant.

## Creation

### Syntax

```matlab
PriorMdl = svar.minnesotainwbvarm(numseries,numlags,Y)
PriorMdl = svar.minnesotainwbvarm(numseries,numlags,Y,Name=Value)
```

### Description

`PriorMdl = svar.minnesotainwbvarm(numseries,numlags,Y)` creates an independent
Normal-Wishart Minnesota prior for a VAR model with `numseries` response
variables and `numlags` autoregressive lags.

`PriorMdl = svar.minnesotainwbvarm(___,Name=Value)` sets residual-variance
handling, Minnesota hyperparameters, and inherited `semiconjugatebvarm` model
options.

## Input Arguments

`numseries` - Number of response series
: Positive integer.

`numlags` - Number of autoregressive lags
: Positive integer.

`Y` - Response data
: Nonempty response sample.

## Name-Value Arguments

`Psi` - Residual-variance scale
: `"exact"` (default) | `"conditional"` | positive numeric vector.

  String values select the estimator used by
  `svar.estimateResidualVariances`. A numeric vector must contain one positive
  value per response series.

`lambda1` - Overall tightness
: `0.2` (default) | positive scalar.

`lambda2` - Cross-variable relative tightness
: `0.5` (default) | positive scalar.

`lambda3` - Lag-decay exponent
: `1` (default) | nonnegative scalar.

`Vc` - Prior variance for deterministic terms
: `1e4` (default) | positive scalar.

`PriorMean` - Own first-lag prior mean
: `ones(1,numseries)` (default) | numeric row vector.

`IncludeConstant` - Flag for including model constant
: Inherited model option.

`IncludeTrend` - Flag for including linear time trend
: Inherited model option.

`NumPredictors` - Number of exogenous predictors
: Inherited model option.

`SeriesNames` - Response series names
: Inherited model option.

`Description` - Model description
: Inherited model option.

## Properties

`ResidualVariances` - Per-series residual variances
: Row vector used as the prior scale.

`lambda1` - Overall Minnesota tightness
: Positive scalar.

`lambda2` - Cross-variable relative tightness
: Positive scalar.

`lambda3` - Lag-decay exponent
: Nonnegative scalar.

`Vc` - Prior variance for deterministic terms
: Positive scalar.

`PriorMean` - Own first-lag prior mean
: Row vector with one value per response series.

Inherited `semiconjugatebvarm` properties such as `NumSeries`, `P`, `Mu`, `V`,
`Omega`, and `DoF` are also available.

## Object Functions

`estimate`
: Estimate the semiconjugate posterior for the stored sample using inherited
  Gibbs sampling. The returned model is an empirical posterior object from
  Econometrics Toolbox.

`simulate`
: Draw reduced-form VAR parameters given the stored sample.

`forecast`
: Forecast responses beyond the stored sample.

Inherited `semiconjugatebvarm` object functions are available where they are
valid for the resulting prior or posterior.

## Examples

### Create an Independent Normal-Wishart Minnesota Prior

```matlab
psi = svar.estimateResidualVariances(Y,4,Method="conditional");
PriorMdl = svar.minnesotainwbvarm(size(Y,2),4,Y,Psi=psi,lambda2=0.4);
```

### Estimate with Reproducible Simulation

Because the semiconjugate posterior is simulation based, set the random number
generator before estimating.

```matlab
rng default
PosteriorMdl = estimate(PriorMdl,Display="off");
```

## More About

### Coefficient Prior Variance

For lag $$\ell$$, source variable $$k$$, target variable $$i$$, and residual
variance scale $$\psi$$, the own-lag target variance is

$$
\operatorname{Var}(B_{\ell,k,i}) =
\frac{\lambda_1^2}{\ell^{2\lambda_3}}\frac{\psi_i}{\psi_k}.
$$

For cross-variable lags, the variance is multiplied by $$\lambda_2^2$$:

$$
\operatorname{Var}(B_{\ell,k,i}) =
\frac{\lambda_1^2\lambda_2^2}{\ell^{2\lambda_3}}\frac{\psi_i}{\psi_k}.
$$

### Scope

This class does not implement sum-of-coefficients or dummy-initial-observation
priors. It also does not provide `logMarginalLikelihood`; use the conjugate
`svar.minnesotamniwbvarm` variant when an analytic marginal likelihood is
required.

## See Also

`minnesotabvarm`, `svar.minnesotamniwbvarm`, `svar.minnesotanbvarm`,
`semiconjugatebvarm`, `svar.estimateResidualVariances`
