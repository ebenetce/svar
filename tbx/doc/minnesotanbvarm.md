# minnesotanbvarm

Fixed-Sigma Normal Minnesota prior for Bayesian VAR models.

The `minnesotanbvarm` model object specifies a Minnesota-shaped Normal prior
with fixed innovations covariance $$\Sigma = \operatorname{diag}(\psi)$$. It
uses the `normalbvarm` parameterization and supports a free cross-variable
tightness hyperparameter, $$\lambda_2$$.

## Creation

### Syntax

```matlab
PriorMdl = minnesotanbvarm(numseries,numlags,ResidualVariances=psi)
PriorMdl = minnesotanbvarm(numseries,numlags,ResidualVariances=psi,Name=Value)
```

### Description

`PriorMdl = minnesotanbvarm(numseries,numlags,ResidualVariances=psi)` creates a
fixed-Sigma Normal Minnesota prior for a VAR model with `numseries` response
variables and `numlags` autoregressive lags.

`PriorMdl = minnesotanbvarm(___,Name=Value)` sets Minnesota hyperparameters and
inherited `normalbvarm` model options.

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

Inherited `normalbvarm` properties such as `NumSeries`, `P`, `Mu`, `V`, and
`Sigma` are also available.

## Object Functions

Inherited `normalbvarm` object functions are available where they are valid for
the resulting prior or posterior.

## Examples

### Create a Fixed-Sigma Minnesota Prior

```matlab
psi = estimateResidualVariances(Y,4,Method="conditional");
PriorMdl = minnesotanbvarm(size(Y,2),4, ...
    ResidualVariances=psi,lambda2=0.4);
```

## More About

### Coefficient Prior Variance

For lag $$\ell$$, source variable $$k$$, target variable $$i$$, and residual
variance scale $$\psi$$, the fixed-Sigma Normal prior uses

$$
\operatorname{Var}(B_{\ell,k,i}) =
\frac{\lambda_1^2}{\ell^{2\lambda_3}}\frac{\psi_i}{\psi_k}
$$

for own lags, and

$$
\operatorname{Var}(B_{\ell,k,i}) =
\frac{\lambda_1^2\lambda_2^2}{\ell^{2\lambda_3}}\frac{\psi_i}{\psi_k}
$$

for cross lags.

## See Also

`minnesotaSpec`, `minnesotanSpec`, `minnesotainwbvarm`, `normalbvarm`,
`estimateResidualVariances`

