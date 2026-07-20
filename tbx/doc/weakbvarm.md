# weakbvarm

Weak improper Bayesian VAR prior for Uhlig-style reduced-form workflows.

The `weakbvarm` model object specifies an improper conjugate Bayesian VAR prior
for reduced-form VAR parameters. It is designed to reproduce the weak
Normal-Wishart prior used in Uhlig (2005), Appendix B, before any structural
sign-restriction or impulse-response identification step is applied.

## Creation

### Syntax

```matlab
PriorMdl = weakbvarm(numseries,numlags)
PriorMdl = weakbvarm(numseries,numlags,Name=Value)
```

### Description

`PriorMdl = weakbvarm(numseries,numlags)` creates a weak improper conjugate
Bayesian VAR prior with `numseries` response variables and `numlags`
autoregressive lags.

`PriorMdl = weakbvarm(numseries,numlags,Name=Value)` sets model options using
name-value arguments. For example, `weakbvarm(6,12,SeriesNames=names)` creates
the six-variable, 12-lag model shape used in Uhlig (2005).

### Input Arguments

`numseries` - Number of response series
: Positive integer.

`numlags` - Number of lagged responses
: Nonnegative integer.

### Name-Value Arguments

`IncludeConstant` - Flag for including model constant
: `true` (default) | `false`

`IncludeTrend` - Flag for including linear time trend
: `false` (default) | `true`

`NumPredictors` - Number of exogenous predictors
: `0` (default) | nonnegative integer

`SeriesNames` - Response series names
: string vector | cell array of character vectors

`Description` - Model description
: string scalar | character vector

## Properties

### Model Characteristics and Dimensionality

`NumSeries` - Number of response series
: Positive integer.

`P` - Multivariate autoregressive polynomial order
: Nonnegative integer.

`SeriesNames` - Response series names
: String vector.

`IncludeConstant` - Flag for including model constant
: Logical scalar.

`IncludeTrend` - Flag for including linear time trend
: Logical scalar.

`NumPredictors` - Number of exogenous predictors
: Nonnegative integer.

`NumParamsEq` - Number of coefficients per equation
: Positive integer. For the default intercept model, this is
  `NumSeries*P + 1`.

### Distribution Hyperparameters

`Mu` - Mean of vectorized matrix normal prior
: Zero vector. The prior is diffuse because `V` has infinite diagonal entries.

`V` - Scaled conditional covariance matrix
: Diagonal matrix with infinite diagonal entries.

`Omega` - Inverse-Wishart scale matrix
: Zero matrix of size `NumSeries`-by-`NumSeries`.

`DoF` - Inverse-Wishart degrees of freedom
: `0`.

## Object Functions

`estimate`
: Estimate the weak-prior reduced-form posterior. The returned posterior is a
  `conjugatebvarm` object.

`simulate`
: Draw reduced-form VAR coefficients and covariance matrices. Use
  `simulate(PriorMdl,Y,NumDraws=N)` to draw from the posterior implied by data
  `Y`, or estimate first and call `simulate(PosteriorMdl,NumDraws=N)`.

Inherited `conjugatebvarm` object functions are also available where they are
valid for the resulting prior or posterior object.

## Examples

### Create Uhlig's Reduced-Form Model Shape

Uhlig (2005) uses monthly data, six variables, 12 lags, and an intercept.

```matlab
PriorMdl = weakbvarm(6,12);
PriorMdl.NumParamsEq
PriorMdl.DoF
```

The number of coefficients per equation is `73`, and the prior degrees of
freedom are `0`.

### Estimate the Reduced-Form Posterior

Estimate a weak-prior posterior from response data `Y`.

```matlab
PriorMdl = weakbvarm(6,12);
PosteriorMdl = estimate(PriorMdl,Y,Display="off");
```

The posterior coefficient mean is the OLS estimate, the coefficient scale is
`inv(X'*X)`, the inverse-Wishart scale is the residual sum of squares, and the
posterior degrees of freedom equal the effective sample size.

### Draw Reduced-Form Parameters

Draw reduced-form VAR parameters directly from the posterior implied by `Y`.

```matlab
PriorMdl = weakbvarm(6,12);
[Coeff,Sigma] = simulate(PriorMdl,Y,NumDraws=1000);
```

You can also estimate first and draw from the posterior object.

```matlab
PosteriorMdl = estimate(PriorMdl,Y,Display="off");
[Coeff,Sigma] = simulate(PosteriorMdl,NumDraws=1000);
```

## More About

### Relation to Uhlig (2005)

Uhlig (2005), Appendix B, describes the reduced-form Normal-Wishart prior using
hyperparameters `B0`, `N0`, `S0`, and `n0`. The weak prior sets

$$
N_0 = 0,\qquad n_0 = 0.
$$

with `S0` and `B0` arbitrary. The resulting posterior is

$$
B_T = \widehat{B},\qquad
S_T = \widehat{S},\qquad
n_T = T,\qquad
N_T = X^{\mathsf{T}}X.
$$

where `Bhat` and `Shat` are the unrestricted reduced-form OLS quantities and
`T` is the effective sample size.

`weakbvarm` encodes the same reduced-form prior in MATLAB's
`conjugatebvarm` parameterization:

$$
\mu = 0,\qquad
V = \operatorname{diag}(\infty),\qquad
\Omega = 0,\qquad
\nu = 0.
$$

### Relation to Mountford and Uhlig (2009)

Mountford and Uhlig (2009) extend Uhlig's sign-restriction approach to fiscal
policy shocks. The paper states that computations use a Bayesian approach as in
Uhlig (2005), drawing VAR coefficients and the covariance matrix from the
posterior before identifying shocks.

This makes Mountford and Uhlig (2009) a useful benchmark for the reduced-form
draw workflow:

- construct the paper's reduced-form VAR data matrix,
- estimate `weakbvarm` with the same lag length and deterministic terms,
- verify posterior moments against direct OLS formulas,
- verify that `simulate` produces coefficient and covariance draws with the
  expected dimensions and positive definite covariance matrices.

It is not, by itself, an independent benchmark for a different reduced-form
prior, because the paper references Uhlig (2005) rather than restating all
Normal-Wishart hyperparameters.

### Scope

`weakbvarm` implements only the reduced-form prior, posterior, and inherited
posterior simulation. It does not identify structural shocks, impose sign or
zero restrictions, compute impulse responses, or construct credible bands.

## See Also

`conjugatebvarm`, `diffusebvarm`, `uniformirbvarm`
