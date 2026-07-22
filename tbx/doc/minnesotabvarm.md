# minnesotabvarm

Create a Minnesota prior for Bayesian VAR models.

`minnesotabvarm` is the public front door for Minnesota priors. It computes
the data-dependent residual-variance scale and dispatches to the requested
Minnesota prior family.

## Syntax

```matlab
PriorMdl = minnesotabvarm(numseries,numlags,Y)
PriorMdl = minnesotabvarm(numseries,numlags,Y,Type=type)
PriorMdl = minnesotabvarm(___,Name=Value)
```

## Description

`PriorMdl = minnesotabvarm(numseries,numlags,Y)` creates a conjugate
Matrix-Normal-Inverse-Wishart Minnesota prior for the sample `Y`.

`PriorMdl = minnesotabvarm(numseries,numlags,Y,Type=type)` selects the
Minnesota prior family. The default is `Type="mniw"`.

`PriorMdl = minnesotabvarm(___,Name=Value)` sets residual-variance handling,
Minnesota hyperparameters, and inherited Bayesian VAR model options.

## Input Arguments

`numseries` - Number of response series
: Positive integer.

`numlags` - Number of autoregressive lags
: Positive integer.

`Y` - Response data
: Nonempty response sample used to estimate residual variances and construct
  data-dependent dummy priors when the selected family uses them.

## Name-Value Arguments

`Type` - Minnesota prior family
: `"mniw"` (default) | `"inw"` | `"normal"`.

  `"mniw"` returns a `svar.minnesotamniwbvarm` object. `"inw"` returns a
  `svar.minnesotainwbvarm` object. `"normal"` returns a
  `svar.minnesotanbvarm` object. Aliases include `"conjugate"` and
  `"matrixnormal"` for `"mniw"`, `"semiconjugate"` and `"kadiyala"` for
  `"inw"`, and `"litterman"` or `"fixedsigma"` for `"normal"`.

`Psi` - Residual-variance scale
: `"exact"` (default) | `"conditional"` | positive numeric vector.

  String values select the estimator used by
  `svar.estimateResidualVariances`. Pass a numeric row vector when tuning
  hyperparameters so the residual-variance scale is not recomputed for every
  candidate.

`lambda1` - Overall tightness
: `0.2` (default) | positive scalar.

`lambda2` - Cross-variable relative tightness
: Positive scalar. Applies to `Type="inw"` and `Type="normal"`.

`lambda3` - Lag-decay exponent
: `1` (default) | nonnegative scalar.

`lambda4` - Sum-of-coefficients dummy tightness
: `Inf` (default) | positive scalar. Applies to `Type="mniw"`.

`lambda5` - Dummy-initial-observation tightness
: `Inf` (default) | positive scalar. Applies to `Type="mniw"`.

`Vc` - Prior variance for deterministic terms
: `1e4` (default) | positive scalar.

`PriorMean` - Own first-lag prior mean
: `ones(1,numseries)` (default) | numeric row vector.

You can also specify inherited Bayesian VAR options such as `SeriesNames`,
`IncludeConstant`, `IncludeTrend`, `NumPredictors`, and `Description`.

## Output Arguments

`PriorMdl` - Minnesota prior model
: `svar.minnesotamniwbvarm`, `svar.minnesotainwbvarm`, or
  `svar.minnesotanbvarm` model object.

## Examples

### Create the Default Conjugate Minnesota Prior

```matlab
numLags = 4;
PriorMdl = minnesotabvarm(size(Y,2),numLags,Y);
PosteriorMdl = estimate(PriorMdl);
```

### Select the Semiconjugate Family

```matlab
PriorMdl = minnesotabvarm(size(Y,2),4,Y,Type="inw",lambda2=0.5);
```

### Reuse Residual Variance Estimates During Tuning

```matlab
psi = svar.estimateResidualVariances(Y,4,Method="conditional");
PriorMdl = minnesotabvarm(size(Y,2),4,Y,Psi=psi,lambda1=0.3);
```

## More About

### Family Dispatch

`minnesotabvarm` is a dispatching function, not a class. Use
`isa(PriorMdl,"svar.minnesotabvarmBase")` to test whether a model is one of the
Minnesota prior families.

## See Also

`svar.minnesotamniwbvarm`, `svar.minnesotainwbvarm`,
`svar.minnesotanbvarm`, `glp`, `svar.estimateResidualVariances`
