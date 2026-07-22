# svar.minnesotamniwbvarm

Conjugate Matrix-Normal-Inverse-Wishart Minnesota prior for Bayesian VAR
models.

The `svar.minnesotamniwbvarm` model object specifies a Litterman-style
Minnesota prior in the `conjugatebvarm` parameterization. The constructor uses
the response sample to resolve residual variances and any active dummy priors.

## Creation

### Syntax

```matlab
PriorMdl = svar.minnesotamniwbvarm(numseries,numlags,Y)
PriorMdl = svar.minnesotamniwbvarm(numseries,numlags,Y,Name=Value)
```

### Description

`PriorMdl = svar.minnesotamniwbvarm(numseries,numlags,Y)` creates a conjugate
Minnesota prior for a VAR model with `numseries` response variables and
`numlags` autoregressive lags. The object stores the sample and the full
dummy-augmented prior.

`PriorMdl = svar.minnesotamniwbvarm(___,Name=Value)` sets residual-variance
handling, Minnesota hyperparameters, and inherited `conjugatebvarm` model
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

`NumDummyObservations` - Number of active dummy observations
: Nonnegative integer.

Inherited `conjugatebvarm` properties such as `NumSeries`, `P`, `SeriesNames`,
`Mu`, `V`, `Omega`, and `DoF` are also available.

## Object Functions

`estimate`
: Estimate the dummy-augmented conjugate posterior for the stored sample. The
  returned posterior is a plain `conjugatebvarm` object.

`simulate`
: Draw reduced-form VAR parameters given the stored sample.

`forecast`
: Forecast responses beyond the stored sample.

`simsmooth`
: Run the inherited simulation smoother using the stored sample.

## Examples

### Create a Minnesota Prior

```matlab
numLags = 4;
psi = svar.estimateResidualVariances(Y,numLags,Method="conditional");
PriorMdl = svar.minnesotamniwbvarm(size(Y,2),numLags,Y,Psi=psi);
```

### Enable Dummy Priors

```matlab
PriorMdl = svar.minnesotamniwbvarm(size(Y,2),4,Y,Psi=psi, ...
    lambda4=10,lambda5=5);
```

### Estimate the Posterior

```matlab
PosteriorMdl = estimate(PriorMdl);
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
data. `svar.minnesotamniwbvarm` folds active dummy priors into the prior at
construction time. Setting `lambda4=Inf` or `lambda5=Inf` disables the
corresponding dummy rows.

## See Also

`minnesotabvarm`, `svar.minnesotainwbvarm`, `svar.minnesotanbvarm`,
`svar.estimateResidualVariances`, `conjugatebvarm`, `glp`
