# glp

Tune conjugate Minnesota hyperparameters by marginal likelihood.

The `glp` function tunes free hyperparameters of
`svar.minnesotamniwbvarm` by maximizing the analytic marginal likelihood for
the response sample. It supports flat bounded searches and embedded
`hyperprior` objects.

## Syntax

```matlab
PriorMdl = glp(numseries,numlags,Y)
PriorMdl = glp(numseries,numlags,Y,Name=Value)
[PriorMdl,info] = glp(___)
```

## Description

`PriorMdl = glp(numseries,numlags,Y)` tunes the default free parameter set for
a conjugate Minnesota prior and returns a `svar.minnesotamniwbvarm` object
built at the optimizer solution.

`PriorMdl = glp(numseries,numlags,Y,Name=Value)` controls fixed and free
Minnesota hyperparameters, residual-variance handling, optimizer options, and
inherited Bayesian VAR model options.

`[PriorMdl,info] = glp(___)` also returns optimization details, including free
parameter names, bounds, starting values, final values, and `fmincon` output.

## Input Arguments

`numseries` - Number of response series
: Positive integer.

`numlags` - Number of autoregressive lags
: Positive integer.

`Y` - Response data
: Numeric matrix, table, or timetable. The number of variables or columns must
  equal `numseries`.

## Name-Value Arguments

`lambda1`, `lambda3`, `lambda4`, `lambda5` - Minnesota hyperparameters
: Scalar | two-element bounds | scalar `hyperprior`.

  A scalar fixes the value. A two-element vector `[lower upper]` makes the
  parameter free inside the bounds with a flat prior. A `hyperprior` makes the
  parameter free inside `hyperprior.Bounds` and adds `hyperprior.logpdf` to the
  objective. `lambda4=Inf` and `lambda5=Inf` disable the corresponding dummy
  priors.

`Vc` - Prior variance for deterministic terms
: `1e4` (default) | positive scalar.

`PriorMean` - Own first-lag prior mean
: `1` (default) | numeric row vector.

`Psi` - Residual-variance scale
: `"exact"` (default) | `"conditional"` | positive numeric vector |
  `2`-by-`numseries` positive bounds | scalar or row vector of `hyperprior`
  objects.

  String values call `svar.estimateResidualVariances` once before the search.
  Numeric vectors fix `Psi`. Numeric bounds or `hyperprior` objects make `Psi`
  free in the optimizer.

`OptimOptions` - Optimization options
: `optimoptions("fmincon",...)` object.

`IncludeConstant` - Flag for including model constant
: Passed to the built `svar.minnesotamniwbvarm` model.

`IncludeTrend` - Flag for including linear time trend
: Passed to the built model.

`NumPredictors` - Number of exogenous predictors
: Passed to the built model.

`SeriesNames` - Response series names
: Passed to the built model. If `Y` is tabular and `SeriesNames` is omitted,
  variable names from `Y` are used.

`Description` - Model description
: Passed to the built model.

## Output Arguments

`PriorMdl` - Tuned Minnesota prior
: `svar.minnesotamniwbvarm` object built from the final hyperparameters and
  residual variances.

`info` - Optimization information
: Structure with fields such as `FreeLambdas`, `PsiNames`, `PsiFree`,
  `FinalLambdas`, `InitialPsi`, `FinalPsi`, `UsedHyperprior`, `X0`,
  `LowerBound`, `UpperBound`, `XHat`, `Objective`, `ExitFlag`, and
  `FminconOutput`.

## Examples

### Tune Minnesota Hyperparameters

```matlab
PriorMdl = glp(size(Y,2),4,Y, ...
    lambda1=hyperprior("Gamma",0.2,0.4,Bounds=[1e-4 5]), ...
    lambda4=hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
    lambda5=hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
    Psi="conditional");
```

### Tune Residual Variances with Hyperpriors

```matlab
psi0 = svar.estimateResidualVariances(Y,4,Method="conditional");
Psi = arrayfun(@(x) hyperprior("InverseGamma",0.02^2,0.02^2, ...
    X0=x,Bounds=[1/100 100]*x), psi0);

[PriorMdl,info] = glp(size(Y,2),4,Y,Psi=Psi);
```

### Tune Within Fixed Bounds

```matlab
PriorMdl = glp(size(Y,2),4,Y, ...
    lambda1=[0.01 1],lambda3=1,lambda4=[1 50],Psi="conditional");
```

## More About

### Objective Function

For free parameters $$\theta$$ and residual variances $$\psi$$, `glp`
minimizes the negative of

$$
\log p(Y \mid \theta,\psi) + \log p(\theta) + \log p(\psi).
$$

The hyperprior terms are included only for fields represented by `hyperprior`
objects. Numeric bounds contribute no density term and therefore behave as flat
priors inside the specified box.

### Supported Prior Family

`glp` tunes `svar.minnesotamniwbvarm`. The analytic marginal likelihood is
available for the conjugate Minnesota family, but not for the independent
Normal-Wishart family.

## See Also

`svar.minnesotamniwbvarm`, `hyperprior`, `logMarginalLikelihood`,
`svar.estimateResidualVariances`, `fmincon`
