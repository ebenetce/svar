# minnesotainwbvarm

Independent Normal-Wishart Minnesota prior for Bayesian VAR models.

The `minnesotainwbvarm` model object specifies a Minnesota prior using the
`semiconjugatebvarm` parameterization. It supports a free cross-variable
tightness hyperparameter, $$\lambda_2$$, at the cost of losing the analytic
marginal likelihood available in the conjugate MNIW variant.

## Creation

### Syntax

```matlab
PriorMdl = minnesotainwbvarm(numseries,numlags,ResidualVariances=psi)
PriorMdl = minnesotainwbvarm(numseries,numlags,ResidualVariances=psi,Name=Value)
```

### Description

`PriorMdl = minnesotainwbvarm(numseries,numlags,ResidualVariances=psi)` creates
an independent Normal-Wishart Minnesota prior for a VAR model with `numseries`
response variables and `numlags` autoregressive lags.

`PriorMdl = minnesotainwbvarm(___,Name=Value)` sets Minnesota hyperparameters
and inherited `semiconjugatebvarm` model options.

## Input Arguments

`numseries` - Number of response series
: Positive integer.

`numlags` - Number of autoregressive lags
: Positive integer.

## Name-Value Arguments

`ResidualVariances` - Per-series residual variances
: Positive numeric vector with `numseries` elements. This argument is required.

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
: Estimate the semiconjugate posterior using inherited Gibbs sampling. The
  returned model is an empirical posterior object from Econometrics Toolbox.

Inherited `semiconjugatebvarm` object functions are available where they are
valid for the resulting prior or posterior.

## Examples

### Create an Independent Normal-Wishart Minnesota Prior

```matlab
psi = estimateResidualVariances(Y,4,Method="conditional");
PriorMdl = minnesotainwbvarm(size(Y,2),4, ...
    ResidualVariances=psi,lambda2=0.4);
```

### Estimate with Reproducible Simulation

Because the semiconjugate posterior is simulation based, set the random number
generator before estimating.

```matlab
rng default
PosteriorMdl = estimate(PriorMdl,Y,Display="off");
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
`minnesotamniwbvarm` variant when an analytic marginal likelihood is required.

## See Also

`minnesotaSpec`, `minnesotamniwbvarm`, `minnesotanbvarm`,
`semiconjugatebvarm`, `estimateResidualVariances`

