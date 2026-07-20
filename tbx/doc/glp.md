# glp

Tune Minnesota hyperparameters by marginal likelihood.

The `glp` function tunes free fields in a conjugate Minnesota specification and
optionally tunes the residual variance scale. It follows the Giannone, Lenza,
and Primiceri style of selecting hyperparameters by maximizing an analytic
marginal likelihood, with optional embedded hyperpriors.

## Syntax

```matlab
mdl = glp(numseries,numlags,Y,Psi)
mdl = glp(numseries,numlags,Y,Psi,Spec=spec)
[mdl,info] = glp(___)
```

## Description

`mdl = glp(numseries,numlags,Y,Psi)` tunes the free hyperparameters in a default
`minnesotaSpec("mniw")` specification and returns a `minnesotamniwbvarm` model.

`mdl = glp(numseries,numlags,Y,Psi,Spec=spec)` uses the supplied
`minnesotamniwSpec` object. Scalar spec fields are fixed, two-element bounds are
tuned with a flat prior, and `hyperprior` fields are tuned with their log-density
added to the objective.

`[mdl,info] = glp(___)` also returns optimization details, including starting
values, bounds, selected hyperparameter names, and the final specification.

## Input Arguments

`numseries` - Number of response series
: Positive integer.

`numlags` - Number of autoregressive lags
: Positive integer.

`Y` - Response data
: Numeric matrix or table. The number of columns must equal `numseries`.

`Psi` - Residual variance specification
: One of these values:

* `1`-by-`numseries` positive numeric vector for fixed residual variances.
* `2`-by-`numseries` numeric matrix of lower and upper bounds.
* Scalar `hyperprior`, broadcast independently across series.
* `1`-by-`numseries` `hyperprior` array.

## Name-Value Arguments

`Spec` - Minnesota specification
: `minnesotaSpec("mniw")` (default) | `minnesotamniwSpec`.

`OptimOptions` - Optimization options
: `optimoptions("fmincon",...)` object. Use this argument to control
  optimizer display, tolerances, and finite-difference settings.

`IncludeConstant` - Flag for including model constant
: Passed to the built `minnesotamniwbvarm` model.

`IncludeTrend` - Flag for including linear time trend
: Passed to the built model.

`NumPredictors` - Number of exogenous predictors
: Passed to the built model.

`SeriesNames` - Response series names
: Passed to the built model. If `Y` is a table and `SeriesNames` is omitted,
  variable names from `Y` are used.

`Description` - Model description
: Passed to the built model.

## Output Arguments

`mdl` - Tuned Minnesota prior
: `minnesotamniwbvarm` object built from the final specification and final
  residual variances.

`info` - Optimization information
: Structure with fields such as `LambdaNames`, `PsiNames`, `InitialSpec`,
  `FinalSpec`, `InitialPsi`, `FinalPsi`, `X0`, `LowerBound`, `UpperBound`,
  `XHat`, `Objective`, `ExitFlag`, and `FminconOutput`.

## Examples

### Tune Minnesota Hyperparameters

Create a specification with embedded hyperpriors and tune both the Minnesota
hyperparameters and residual variances.

```matlab
spec = minnesotaSpec("mniw", ...
    lambda1=hyperprior("Gamma",0.2,0.4,Bounds=[1e-4 5]), ...
    lambda4=hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
    lambda5=hyperprior("Gamma",1,1,Bounds=[1e-4 50]));

psi0 = estimateResidualVariances(Y,1,Method="conditional");
Psi = arrayfun(@(x) hyperprior("InverseGamma",0.02^2,0.02^2, ...
    X0=x, Bounds=[1/100 100]*x), psi0);

[PriorMdl,info] = glp(size(Y,2),4,Y,Psi,Spec=spec);
```

### Tune Within Fixed Bounds

Use numeric bounds for a flat-prior search over selected fields.

```matlab
spec = minnesotaSpec("mniw",lambda1=[0.01 1],lambda3=1, ...
    lambda4=[1 50],lambda5=Inf);

psi0 = estimateResidualVariances(Y,4,Method="conditional");
PriorMdl = glp(size(Y,2),4,Y,psi0,Spec=spec);
```

## More About

### Objective Function

For free parameters $$\theta$$ and residual variances $$\psi$$, `glp` minimizes
the negative of the marginal-likelihood objective

$$
\log p(Y \mid \theta,\psi) + \log p(\theta) + \log p(\psi).
$$

The hyperprior terms are included only for fields represented by `hyperprior`
objects. Numeric bounds contribute no density term and therefore behave as flat
priors inside the specified box.

### Supported Specification Family

`glp` requires `minnesotaSpec("mniw")`. The analytic marginal likelihood is
available for the conjugate MNIW Minnesota prior, but not for the independent
Normal-Wishart or fixed-Sigma Normal variants.

## See Also

`minnesotaSpec`, `minnesotamniwbvarm`, `hyperprior`, `estimateResidualVariances`,
`fmincon`

