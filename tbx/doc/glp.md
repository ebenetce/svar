# glp

Tune Minnesota hyperparameters by marginal likelihood, and optionally sample them.

The `glp` function chooses the free hyperparameters of a conjugate Minnesota
prior by maximizing its analytic marginal likelihood, optionally combined with
hyperpriors. It implements the approach of Giannone, Lenza, and Primiceri
(2012), including the Metropolis step that integrates over the hyperparameters
rather than fixing them at the maximizer.

## Syntax

```matlab
mdl = glp(numseries,numlags,Y)
mdl = glp(numseries,numlags,Y,Name=Value)
[mdl,info] = glp(___)
[mdl,info,chain] = glp(___,NumDraws=m)
```

## Description

`mdl = glp(numseries,numlags,Y)` tunes the free hyperparameters against `Y` and
returns a `svar.minnesotamniwbvarm` prior built at the maximizer.

Each hyperparameter is in one of three states, so the free/fixed split is
implicit in what you pass and there is no separate list to keep in sync:

| Value | State |
| --- | --- |
| scalar | Fixed at that value |
| `[lower upper]` | Free, tuned within those bounds under a flat prior |
| `hyperprior` | Free, tuned within its `Bounds` under its log density |

`[mdl,info] = glp(___)` also returns a structure describing the search.

`[mdl,info,chain] = glp(___,NumDraws=m)` additionally samples the hyperparameter
posterior and returns `m` draws of the hyperparameters together with the VAR
coefficients and innovations covariance they imply.

## Input Arguments

`numseries` - Number of response series
: Positive integer.

`numlags` - Number of autoregressive lags
: Positive integer.

`Y` - Response data
: Numeric matrix, table, or timetable with `numseries` columns.

## Name-Value Arguments

`lambda1` - Overall Minnesota tightness
: `0.2` (default) | scalar | `[lower upper]` | `hyperprior`. The paper's
  $$\lambda$$.

`lambda3` - Lag-decay exponent
: `1` (default) | scalar | `[lower upper]` | `hyperprior`. Lag variance decays
  as $$1/l^{2\lambda_3}$$, so `lambda3` is half the paper's $$\alpha$$.

`lambda4` - Sum-of-coefficients tightness
: `Inf` (default, prior off) | scalar | `[lower upper]` | `hyperprior`. The
  paper's $$\mu$$.

`lambda5` - Dummy-initial-observation tightness
: `Inf` (default, prior off) | scalar | `[lower upper]` | `hyperprior`. The
  paper's $$\delta$$.

`Vc` - Prior variance of the constant and trend
: `1e4` (default) | positive scalar. Always fixed.

`PriorMean` - Minnesota prior mean on the own first lag
: `1` (default) | `1`-by-`numseries` vector. Always fixed. Use `0` for series
  entering in first differences.

`Psi` - Residual variance scale
: `"exact"` (default) | `"conditional"` | numeric | `hyperprior`. One of:

* `"exact"` or `"conditional"` - estimated from `Y` once, then held fixed. The
  string form is resolved before the search, never per candidate.
* `1`-by-`numseries` numeric - fixed at these variances.
* `2`-by-`numseries` numeric - free within these `[lower; upper]` bounds.
* scalar `hyperprior` - free, broadcast independently to all series.
* `1`-by-`numseries` `hyperprior` array - free, one log density per series.

`OptimOptions` - Optimization options
: `optimoptions("fmincon",...)` object controlling display, tolerances, and
  finite-difference settings.

`NumDraws` - Hyperparameter draws to keep
: `0` (default) | nonnegative integer. `0` maximizes only. Any positive value
  runs the Metropolis sampler and populates `chain`.

`BurnIn` - Draws discarded before keeping any
: `NumDraws` (default) | nonnegative integer.

`ProposalScale` - Random-walk step size
: `1` (default) | positive scalar. The proposal covariance is
  `ProposalScale^2` times the inverse Hessian at the maximizer. Tune for an
  acceptance rate of roughly 0.2 to 0.3; $$2.38/\sqrt{d}$$ for `d` free
  hyperparameters is a good starting point.

`IncludeConstant`, `IncludeTrend`, `NumPredictors`, `SeriesNames`, `Description`
: Passed to the built `svar.minnesotamniwbvarm` model. If `Y` is tabular and
  `SeriesNames` is omitted, its variable names are used.

## Output Arguments

`mdl` - Tuned Minnesota prior
: `svar.minnesotamniwbvarm` object built at the maximizer.

`info` - Search information
: Structure with fields `FreeLambdas`, `PsiNames`, `PsiFree`, `FinalLambdas`,
  `InitialPsi`, `FinalPsi`, `UsedHyperprior`, `X0`, `LowerBound`, `UpperBound`,
  `XHat`, `Objective`, `ExitFlag`, `FminconOutput`, `NumDraws`,
  `ProposalScale`, and `Hessian`.

`chain` - Hyperparameter posterior draws
: Empty structure when `NumDraws` is `0`. Otherwise a structure with fields:

| Field | Size | Contents |
| --- | --- | --- |
| `lambda1`, `lambda3`, `lambda4`, `lambda5` | `NumDraws`-by-`1` | Hyperparameter draws; fixed ones are constant |
| `Psi` | `NumDraws`-by-`numseries` | Residual variance draws |
| `Coefficients` | `m`-by-`numseries`-by-`NumDraws` | VAR coefficient draws |
| `Sigma` | `numseries`-by-`numseries`-by-`NumDraws` | Innovations covariance draws |
| `LogPosterior` | `NumDraws`-by-`1` | Log posterior at each draw |
| `AcceptanceRate` | scalar | Fraction of proposals accepted |
| `ProposalCovariance` | `d`-by-`d` | Covariance actually used |

## Examples

### Tune Minnesota Hyperparameters

```matlab
psi0 = estimateResidualVariances(Y,1,Method="conditional");
Psi = arrayfun(@(x) hyperprior("InverseGamma",0.02^2,0.02^2, ...
    Parameterization="native", X0=x, Bounds=[x/100 x*100]), psi0);

[PriorMdl,info] = glp(size(Y,2),5,Y, ...
    lambda1=hyperprior("Gamma",0.2,0.4,Bounds=[1e-4 5]), ...
    lambda4=hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
    lambda5=hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
    Vc=10e6, Psi=Psi);
```

### Integrate Over the Hyperparameters

Bands built from `chain` account for uncertainty about the hyperparameters;
bands built from `PriorMdl` alone condition on a single value of them.

```matlab
[PriorMdl,info,chain] = glp(size(Y,2),5,Y, ...
    lambda1=hyperprior("Gamma",0.2,0.4,Bounds=[1e-4 5]), ...
    Psi=Psi, NumDraws=10000, BurnIn=10000, ProposalScale=0.8);

fprintf("acceptance rate %.3f\n",chain.AcceptanceRate);
histogram(chain.lambda1);
```

### Tune Within Fixed Bounds

Numeric bounds give a flat-prior search over the specified box.

```matlab
psi0 = estimateResidualVariances(Y,5,Method="conditional");
PriorMdl = glp(size(Y,2),5,Y,lambda1=[0.01 1],lambda3=1, ...
    lambda4=[1 50],lambda5=Inf,Psi=psi0);
```

## More About

### Objective Function

For free hyperparameters $$\theta$$ and residual variances $$\psi$$, `glp`
minimizes the negative of

$$
\log p(Y \mid \theta,\psi) + \log p(\theta) + \log p(\psi).
$$

Hyperprior terms enter only for fields represented by `hyperprior` objects.
Numeric bounds contribute no density term and so act as flat priors inside the
specified box.

### Sampling the Hyperparameters

With `NumDraws`, `glp` runs the Metropolis algorithm of the paper's appendix B.
The chain starts at the maximizer, proposes from a Gaussian centered on the
current draw whose covariance is the scaled inverse Hessian, and accepts on the
log posterior ratio. The proposal is symmetric, so no Hastings correction
applies, and candidates outside the search box are rejected rather than clipped.

The Hessian is computed by central differences on the objective itself, not
taken from `fmincon`. The optimizer's quasi-Newton approximation is an estimate
of the Hessian of the Lagrangian, contaminated by barrier terms, and is
generally too ill-conditioned to invert into a usable proposal covariance.

Because the search runs in natural, box-constrained coordinates, no Jacobian
correction is needed - unlike implementations that maximize in a transformed
unconstrained space and must map the Hessian back.

Conditional on each hyperparameter draw, the VAR coefficients and innovations
covariance are drawn from their exact Normal-Inverse-Wishart posterior, so the
kept draws are independent given the hyperparameters.

### Supported Prior Family

`glp` builds `svar.minnesotamniwbvarm`. The analytic marginal likelihood exists
for the conjugate MNIW Minnesota prior, but not for the independent
Normal-Wishart or fixed-Sigma Normal variants.

## See Also

`svar.minnesotamniwbvarm`, `hyperprior`, `logMarginalLikelihood`,
`estimateResidualVariances`, `fmincon`
