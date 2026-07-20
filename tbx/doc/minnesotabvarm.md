# minnesotabvarm

Abstract shared base class for Minnesota BVAR priors.

`minnesotabvarm` stores the hyperparameters and shared helper behavior used by
the concrete Minnesota prior classes. This class is abstract and hidden;
construct `minnesotamniwbvarm`, `minnesotainwbvarm`, or `minnesotanbvarm`
instead.

## Class Details

`minnesotabvarm` is mixed into the concrete Minnesota model classes. It is not
an Econometrics Toolbox model superclass by itself. The concrete classes pair it
with the appropriate Bayesian VAR superclass:

* `minnesotamniwbvarm` uses `conjugatebvarm`.
* `minnesotainwbvarm` uses `semiconjugatebvarm`.
* `minnesotanbvarm` uses `normalbvarm`.

## Properties

`ResidualVariances` - Per-series residual variances
: Positive row vector used as the Minnesota variance scale.

`lambda1` - Overall tightness
: Positive scalar.

`lambda3` - Lag-decay exponent
: Nonnegative scalar.

`Vc` - Prior variance for deterministic terms
: Positive scalar.

`PriorMean` - Own first-lag prior mean
: Row vector with one value per response series.

## Dependent Properties

`m` - Number of coefficients per equation
: Total regressors per equation, equal to lagged response coefficients plus
  deterministic and predictor regressors.

`nex` - Number of deterministic and predictor regressors
: Sum of the constant, trend, and exogenous predictor counts.

## More About

### Shared Prior Mean

The base class builds the Minnesota prior mean with nonzero entries only on own
first lags. For target series $$i$$, the own first-lag prior mean is
`PriorMean(i)`.

### Input Validation

Concrete subclasses use `validateMinnesotaInputs` to require one residual
variance per response series and to resolve an empty `PriorMean` to a row vector
of ones.

## See Also

`minnesotamniwbvarm`, `minnesotainwbvarm`, `minnesotanbvarm`,
`minnesotaSpec`, `estimateResidualVariances`

