# minnesotanSpec

Hyperparameter recipe for a fixed-Sigma Normal Minnesota prior.

`minnesotanSpec` stores the Minnesota hyperparameters used to build a
`minnesotanbvarm` model. Construct this class with `minnesotaSpec("normal")`.

## Creation

### Syntax

```matlab
spec = minnesotaSpec("normal")
spec = minnesotaSpec("normal",Name=Value)
```

### Description

`spec = minnesotaSpec("normal")` creates a fixed-Sigma Normal Minnesota
specification with default hyperparameters.

`spec = minnesotaSpec("normal",Name=Value)` sets hyperparameter fields. Scalar
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
: Materialize a `minnesotanbvarm` model from a resolved specification and a
  residual variance vector.

`pack`
: Return optimizer starts, lower bounds, upper bounds, and field names for
  free hyperparameters.

`unpack`
: Write optimizer values back to selected hyperparameter fields.

`freeFields`
: Return names of free hyperparameters.

## Examples

### Build a Fixed-Sigma Minnesota Prior

```matlab
psi = estimateResidualVariances(Y,4,Method="conditional");
spec = minnesotaSpec("normal",lambda1=0.2,lambda2=0.5,lambda3=1);

PriorMdl = spec.build(size(Y,2),4,psi);
```

## More About

### Fixed Innovations Covariance

The model built from this spec fixes the innovations covariance matrix to

$$
\Sigma = \operatorname{diag}(\psi).
$$

The coefficient covariance is diagonal over coefficient-target pairs, so a free
$$\lambda_2$$ can control cross-variable shrinkage.

## See Also

`minnesotaSpec`, `minnesotanbvarm`, `minnesotainwSpec`,
`estimateResidualVariances`

