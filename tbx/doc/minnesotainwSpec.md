# minnesotainwSpec

Hyperparameter recipe for an independent Normal-Wishart Minnesota prior.

`minnesotainwSpec` stores the Minnesota hyperparameters used to build a
`minnesotainwbvarm` model. Construct this class with `minnesotaSpec("inw")`.

## Creation

### Syntax

```matlab
spec = minnesotaSpec("inw")
spec = minnesotaSpec("inw",Name=Value)
```

### Description

`spec = minnesotaSpec("inw")` creates an independent Normal-Wishart Minnesota
specification with default hyperparameters.

`spec = minnesotaSpec("inw",Name=Value)` sets hyperparameter fields. Scalar
fields are fixed; two-element bounds are free fields for an external search.

## Name-Value Arguments

`lambda1` - Overall tightness
: `0.2` (default) | positive scalar | `[lower upper]` bounds |
  `hyperprior`.

`lambda2` - Cross-variable relative tightness
: `0.5` (default) | positive scalar | `[lower upper]` bounds |
  `hyperprior`.

`lambda3` - Lag-decay exponent
: `1` (default) | nonnegative scalar | `[lower upper]` bounds |
  `hyperprior`.

`Vc` - Prior variance for deterministic terms
: `1e4` (default) | positive scalar.

`PriorMean` - Own first-lag prior mean
: `[]` (default) | numeric row vector.

## Properties

`lambda1` - Overall tightness
: Positive scalar, bounds, or `hyperprior`.

`lambda2` - Cross-variable relative tightness
: Positive scalar, bounds, or `hyperprior`.

`lambda3` - Lag-decay exponent
: Nonnegative scalar, bounds, or `hyperprior`.

`Vc` - Prior variance for deterministic terms
: Positive scalar.

`PriorMean` - Own first-lag prior mean
: Numeric row vector or `[]`.

## Object Functions

`build`
: Materialize a `minnesotainwbvarm` model from a resolved specification and a
  residual variance vector.

`pack`
: Return optimizer starts, lower bounds, upper bounds, and field names for
  free hyperparameters.

`unpack`
: Write optimizer values back to selected hyperparameter fields.

`freeFields`
: Return names of free hyperparameters.

## Examples

### Create and Build an INW Prior

```matlab
psi = estimateResidualVariances(Y,4,Method="conditional");
spec = minnesotaSpec("inw",lambda1=0.2,lambda2=0.5,lambda3=1);

PriorMdl = spec.build(size(Y,2),4,psi, ...
    SeriesNames=["Output" "Prices" "Rate"]);
```

### Pack Free Fields for an External Search

```matlab
spec = minnesotaSpec("inw",lambda1=[0.01 1],lambda2=[0.05 1]);
[x0,lb,ub,names] = spec.pack();
candidate = spec.unpack(x0,names);
```

## More About

### Independent Normal-Wishart Prior

The INW variant has a full diagonal coefficient covariance over
coefficient-target pairs. This allows separate own- and cross-variable
shrinkage:

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

### Tuning

This specification does not have the analytic marginal likelihood used by
`glp`. Use an out-of-sample criterion, grid search, or another external
selection method when tuning `lambda2` or other INW fields.

## See Also

`minnesotaSpec`, `minnesotainwbvarm`, `minnesotaBaseSpec`,
`estimateResidualVariances`

