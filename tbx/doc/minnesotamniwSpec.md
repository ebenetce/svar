# minnesotamniwSpec

Hyperparameter recipe for a conjugate MNIW Minnesota prior.

`minnesotamniwSpec` stores the Minnesota hyperparameters used to build a
`minnesotamniwbvarm` model. Construct this class with `minnesotaSpec` or
`minnesotaSpec("mniw")`.

## Creation

### Syntax

```matlab
spec = minnesotaSpec
spec = minnesotaSpec("mniw")
spec = minnesotaSpec("mniw",Name=Value)
```

### Description

`spec = minnesotaSpec` creates the default conjugate MNIW Minnesota
specification.

`spec = minnesotaSpec("mniw",Name=Value)` sets hyperparameter fields. Scalar
fields are fixed, two-element bounds are free with flat bounds, and `hyperprior`
fields are free with a log-prior contribution.

## Name-Value Arguments

`lambda1` - Overall tightness
: `0.2` (default) | positive scalar | `[lower upper]` bounds |
  `hyperprior`.

`lambda3` - Lag-decay exponent
: `1` (default) | nonnegative scalar | `[lower upper]` bounds |
  `hyperprior`.

`lambda4` - Sum-of-coefficients dummy tightness
: `Inf` (default) | positive scalar | `[lower upper]` bounds |
  `hyperprior`. Set to `Inf` to disable this dummy prior.

`lambda5` - Dummy-initial-observation tightness
: `Inf` (default) | positive scalar | `[lower upper]` bounds |
  `hyperprior`. Set to `Inf` to disable this dummy prior.

`Vc` - Prior variance for deterministic terms
: `1e4` (default) | positive scalar.

`PriorMean` - Own first-lag prior mean
: `[]` (default) | numeric row vector.

## Properties

`lambda1` - Overall tightness
: Positive scalar, bounds, or `hyperprior`.

`lambda3` - Lag-decay exponent
: Nonnegative scalar, bounds, or `hyperprior`.

`lambda4` - Sum-of-coefficients dummy tightness
: Positive scalar, `Inf`, bounds, or `hyperprior`.

`lambda5` - Dummy-initial-observation tightness
: Positive scalar, `Inf`, bounds, or `hyperprior`.

`Vc` - Prior variance for deterministic terms
: Positive scalar.

`PriorMean` - Own first-lag prior mean
: Numeric row vector or `[]`.

## Object Functions

`build`
: Materialize a `minnesotamniwbvarm` model from a resolved specification and a
  residual variance vector.

`pack`
: Return optimizer starts, lower bounds, upper bounds, and field names for
  free hyperparameters.

`unpack`
: Write optimizer values back to selected hyperparameter fields.

`logHyperprior`
: Evaluate log-density terms from embedded `hyperprior` fields.

`freeFields`
: Return names of free hyperparameters.

## Examples

### Build a Conjugate Minnesota Prior

```matlab
psi = estimateResidualVariances(Y,4,Method="conditional");
spec = minnesotaSpec("mniw",lambda1=0.25,lambda4=10,lambda5=5);

PriorMdl = spec.build(size(Y,2),4,psi, ...
    SeriesNames=["Output" "Prices" "Rate"]);
```

### Pack Free Fields

```matlab
spec = minnesotaSpec("mniw",lambda1=[0.01 1],lambda4=[1 50]);
[x0,lb,ub,names] = spec.pack();
candidate = spec.unpack(x0,names);
```

## More About

### Conjugate Structure

The MNIW prior preserves conjugacy by using a Kronecker covariance structure.
The implied lag-coefficient variance is

$$
\operatorname{Var}(B_{\ell,k,i}) =
\frac{\lambda_1^2}{\ell^{2\lambda_3}}\frac{\psi_i}{\psi_k}.
$$

This structure pins the cross-variable tightness parameter to
$$\lambda_2 = 1$$. Use `minnesotainwSpec` or `minnesotanSpec` when a free
$$\lambda_2$$ is required.

### Dummy Priors

Finite values of $$\lambda_4$$ and $$\lambda_5$$ activate the
sum-of-coefficients and dummy-initial-observation priors. `Inf` disables the
corresponding dummy rows.

## See Also

`minnesotaSpec`, `minnesotamniwbvarm`, `minnesotainwSpec`, `minnesotanSpec`,
`hyperprior`, `glp`

