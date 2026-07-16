# uniformirbvarm

Bayesian VAR prior with the reduced-form prior induced by uniform impulse
response parameters.

The `uniformirbvarm` model object specifies an improper conjugate Bayesian VAR
prior for the reduced-form parameters of an `n`-D VAR(`p`) model. The class is a
small wrapper around `conjugatebvarm`; it sets the coefficient prior to diffuse,
the inverse-Wishart scale matrix to zero, and the inverse-Wishart degrees of
freedom to the value implied by the uniform impulse-response prior in Arias,
Rubio-Ramirez, and Waggoner.

## Creation

### Syntax

```matlab
PriorMdl = uniformirbvarm(numseries,numlags)
PriorMdl = uniformirbvarm(numseries,numlags,Name=Value)
```

### Description

`PriorMdl = uniformirbvarm(numseries,numlags)` creates a reduced-form Bayesian
VAR prior with `numseries` response variables and `numlags` autoregressive lags.

`PriorMdl = uniformirbvarm(numseries,numlags,Name=Value)` sets model options
using name-value arguments. For example,
`uniformirbvarm(4,4,DeterminantShift=-3)` uses the paper's default determinant
shift.

### Input Arguments

`numseries` - Number of response series
: Positive integer.

`numlags` - Number of lagged responses
: Nonnegative integer.

### Name-Value Arguments

`DeterminantShift` - Determinant exponent shift
: `-3` (default) | numeric scalar

  Shift applied to the number of coefficients per equation in the reduced-form
  determinant exponent. If `k = NumEquationCoefficients`, then the prior density
  is proportional to

  ```text
  |det(Sigma)|^((k + DeterminantShift)/2).
  ```

  Set `DeterminantShift=-3` to use the paper's uniform impulse-response prior.
  This property can be set only when you create the object.

`IncludeConstant` - Flag for including model constant
: `true` (default) | `false`

`IncludeTrend` - Flag for including linear time trend
: `false` (default) | `true`

`NumPredictors` - Number of exogenous predictors
: `0` (default) | nonnegative integer

`SeriesNames` - Response series names
: string vector | cell array of character vectors

Other `conjugatebvarm` name-value arguments can be passed to the superclass, but
`uniformirbvarm` overwrites the prior hyperparameters required by the reduced-form
uniform impulse-response prior.

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

`NumEquationCoefficients` - Number of coefficients per equation
: Positive integer. This is `NumSeries*P + IncludeConstant + IncludeTrend +
  NumPredictors`.

### Uniform Prior Parameters

`DeterminantShift` - Determinant exponent shift
: Numeric scalar. The default is `-3`. Set access is private.

`ReducedFormLogDetExponent` - Reduced-form log-determinant exponent
: Numeric scalar equal to `(NumEquationCoefficients + DeterminantShift)/2`.

### Distribution Hyperparameters

`Mu` - Mean of vectorized matrix normal prior
: Zero vector of length `NumSeries*NumEquationCoefficients`.

`V` - Scaled conditional covariance matrix
: Diagonal matrix with infinite diagonal entries.

`Omega` - Inverse-Wishart scale matrix
: Zero matrix of size `NumSeries`-by-`NumSeries`.

`DoF` - Inverse-Wishart degrees of freedom
: Numeric scalar equal to
  `-2*NumEquationCoefficients - NumSeries - 1 - DeterminantShift`.

## Object Functions

`estimate`
: Estimate the conjugate reduced-form posterior. The returned posterior is a
  `conjugatebvarm` object.

Inherited `conjugatebvarm` object functions are also available where they are
valid for the resulting prior or posterior object.

## Examples

### Create the Paper Reduced-Form Prior

Create the reduced-form prior for a four-variable VAR with four lags and an
intercept.

```matlab
PriorMdl = uniformirbvarm(4,4)
```

The model uses `DeterminantShift=-3`. Since `k = 4*4 + 1 = 17`, the
reduced-form determinant exponent is `7` and the prior degrees of freedom are
`-36`.

```matlab
PriorMdl.NumEquationCoefficients
PriorMdl.ReducedFormLogDetExponent
PriorMdl.DoF
```

### Estimate a Reduced-Form Posterior

Estimate the posterior from a matrix of responses `Y`.

```matlab
PriorMdl = uniformirbvarm(4,4,SeriesNames=["OPHNFB" "Hours" "GDPDEF" "GS10"]);
PosteriorMdl = estimate(PriorMdl,Y,Display="off");
```

### Use a Different Determinant Shift

Change the determinant shift when constructing the prior.

```matlab
PriorMdl = uniformirbvarm(4,4,DeterminantShift=-1);
```

Changing the determinant shift changes both `ReducedFormLogDetExponent` and
`DoF`.

## More About

### Relation to the Paper

Proposition 5 in Arias, Rubio-Ramirez, and Waggoner writes a generic
reduced-form prior density as

$$\left|det(\Sigma)\right|^{\frac{a}{2}}$$

Corollary 1 specializes the generic numerator to `a = m - 3`, where `m` is the
number of coefficients per equation. In `uniformirbvarm`, this is represented as

```text
a = NumEquationCoefficients + DeterminantShift
```

with `DeterminantShift=-3` by default.

### Scope

`uniformirbvarm` implements only the reduced-form prior and posterior algebra.
It does not draw orthogonal matrices, impose sign or zero restrictions, compute
structural impulse responses, or construct joint credible sets.

## See Also

`conjugatebvarm`, `diffusebvarm`, `weakbvarm`
