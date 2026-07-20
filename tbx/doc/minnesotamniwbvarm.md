# minnesotamniwbvarm

Conjugate Matrix-Normal-Inverse-Wishart Minnesota prior for Bayesian VAR models.

The `minnesotamniwbvarm` model object specifies a Litterman-style Minnesota
prior in the `conjugatebvarm` parameterization. The constructor is data-free:
it requires precomputed per-series residual variances and stores only the prior
hyperparameters and reduced-form model layout.

## Creation

### Syntax

```matlab
PriorMdl = minnesotamniwbvarm(numseries,numlags,ResidualVariances=psi)
PriorMdl = minnesotamniwbvarm(numseries,numlags,ResidualVariances=psi,Name=Value)
```

### Description

`PriorMdl = minnesotamniwbvarm(numseries,numlags,ResidualVariances=psi)` creates
a conjugate Minnesota prior for a VAR model with `numseries` response variables
and `numlags` autoregressive lags. `psi` is a row vector of residual variance
estimates, one per response series.

`PriorMdl = minnesotamniwbvarm(___,Name=Value)` sets Minnesota hyperparameters
and inherited `conjugatebvarm` model options.

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

`lambda3` - Lag-decay exponent
: `1` (default) | nonnegative scalar.

`lambda4` - Sum-of-coefficients dummy tightness
: `Inf` (default) | positive scalar. Set to `Inf` to disable this dummy prior.

`lambda5` - Dummy-initial-observation tightness
: `Inf` (default) | positive scalar. Set to `Inf` to disable this dummy prior.

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
: Fixed at `1`.

`lambda3` - Lag-decay exponent
: Nonnegative scalar.

`lambda4` - Sum-of-coefficients dummy tightness
: Positive scalar or `Inf`.

`lambda5` - Dummy-initial-observation tightness
: Positive scalar or `Inf`.

`Vc` - Prior variance for deterministic terms
: Positive scalar.

`PriorMean` - Own first-lag prior mean
: Row vector with one value per response series.

Inherited `conjugatebvarm` properties such as `NumSeries`, `P`, `SeriesNames`,
`Mu`, `V`, `Omega`, and `DoF` are also available.

## Object Functions

`estimate`
: Estimate the dummy-augmented conjugate posterior. The returned posterior is a
  plain `conjugatebvarm` object.

`logMarginalLikelihood`
: Evaluate the analytic log marginal likelihood for the hyperparameters.

`negativeLogMarginalLikelihood`
: Return the negative log marginal likelihood for use as an optimizer
  objective.

`marginalLikelihood`
: Return the marginal likelihood on the original scale.

`simulate`
: Draw reduced-form VAR parameters from the posterior implied by data and the
  dummy-augmented prior.

`forecast`
: Forecast from the posterior implied by data and the dummy-augmented prior.

`simsmooth`
: Run the inherited simulation smoother using the dummy-augmented prior.

## Examples

### Create a Minnesota Prior

```matlab
numLags = 4;
psi = estimateResidualVariances(Y,numLags,Method="conditional");
PriorMdl = minnesotamniwbvarm(size(Y,2),numLags, ...
    ResidualVariances=psi);
```

### Enable Dummy Priors

```matlab
PriorMdl = minnesotamniwbvarm(size(Y,2),4,ResidualVariances=psi, ...
    lambda4=10,lambda5=5);
```

### Evaluate a Hyperparameter Objective

```matlab
obj = PriorMdl.negativeLogMarginalLikelihood(Y);
```

## More About

### Conjugate Minnesota Variance

For lag $$\ell$$, source variable $$k$$, target variable $$i$$, and residual
variance scale $$\psi$$, the implied coefficient variance is

$$
\operatorname{Var}(B_{\ell,k,i}) =
\frac{\lambda_1^2}{\ell^{2\lambda_3}}\frac{\psi_i}{\psi_k}.
$$

The conjugate Kronecker structure keeps $$\lambda_2$$ fixed at `1`. This is the
tradeoff that makes the analytic marginal likelihood available.

### Dummy Priors

The sum-of-coefficients and dummy-initial-observation priors depend on the
data. `minnesotamniwbvarm` applies them inside `estimate`, `simulate`,
`forecast`, and the marginal-likelihood methods. Setting `lambda4=Inf` or
`lambda5=Inf` disables the corresponding dummy rows.

## See Also

`minnesotaSpec`, `minnesotamniwSpec`, `minnesotainwbvarm`, `minnesotanbvarm`,
`estimateResidualVariances`, `conjugatebvarm`, `glp`

