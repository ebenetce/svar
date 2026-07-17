# minnesotaSpec

Factory for Minnesota BVAR hyperparameter specifications.

The `minnesotaSpec` dispatcher creates Minnesota-prior spec objects. With no
method argument it returns a `minnesotamniwSpec`, preserving the GLP-style
conjugate workflow.

## Creation

### Syntax

```matlab
spec = minnesotaSpec
spec = minnesotaSpec(Method)
spec = minnesotaSpec(Name=Value)
```

### Description

`spec = minnesotaSpec` creates a `minnesotamniwSpec` with default Minnesota
hyperparameters.

`spec = minnesotaSpec(Method)` chooses `"mniw"`, `"inw"`, or `"normal"`.

`spec = minnesotaSpec(Name=Value)` sets one or more hyperparameter fields.

### Name-Value Arguments

`lambda1` - Overall tightness
: `0.2` (default) | positive scalar.

`lambda3` - Lag-decay exponent
: `1` (default) | nonnegative scalar.

`lambda4` - Sum-of-coefficients dummy tightness
: `Inf` (default) | positive scalar. Set to `Inf` to disable this dummy prior.

`lambda5` - Dummy-initial-observation tightness
: `Inf` (default) | positive scalar. Set to `Inf` to disable this dummy prior.

`Vc` - Prior variance for deterministic terms
: `1e4` (default) | positive scalar.

`PriorMean` - Own first-lag prior mean
: `[]` (default) | numeric row vector. An empty value means
  the built model uses a row vector of ones.

## Properties

`lambda1` - Overall Minnesota tightness
: Positive scalar.

`lambda3` - Lag-decay exponent
: Nonnegative scalar.

`lambda4` - Sum-of-coefficients dummy tightness
: Positive scalar or `Inf`.

`lambda5` - Dummy-initial-observation tightness
: Positive scalar or `Inf`.

`Vc` - Prior variance for constants, trends, and predictors
: Positive scalar.

`PriorMean` - Own first-lag prior mean
: Numeric row vector or `[]`.

## Object Functions

`build`
: Materialize a concrete Minnesota model from the spec and a residual variance
  vector.

`pack`
: Map selected finite positive hyperparameters to log space for unconstrained
  optimization.

`unpack`
: Write an unconstrained optimizer vector back to selected hyperparameter
  fields, returning a modified copy of the spec.

`logHyperprior`
: Evaluate optional Gamma and inverse-Gamma log-density terms for
  hyperparameter tuning.

## Examples

### Build a Prior from a Spec

```matlab
numLags = 4;
psi = estimateResidualVariances(Y,numLags,Method="conditional");

spec = minnesotaSpec(lambda1=0.25,lambda4=10,lambda5=5);
PriorMdl = spec.build(size(Y,2),numLags,psi, ...
    SeriesNames=["Output" "Prices" "Rate"]);
```

### Pack and Unpack Free Parameters

```matlab
freeParams = ["lambda1" "lambda4" "lambda5"];
theta0 = spec.pack(freeParams);

candidate = spec.unpack(theta0 + [0.1 0 0],freeParams);
```

### Add a Hyperprior Term

```matlab
priorcoef.lambda1.k = 2;
priorcoef.lambda1.theta = 0.2;

logp = spec.logHyperprior(priorcoef,psi);
```

## More About

### Separation from Model Objects

`minnesotaSpec` is the mutable recipe. It stores hyperparameter values and
optimization helpers but no data and no model moments. Concrete model objects
such as `minnesotamniwbvarm` are materialized from a spec, a model size, a lag
order, and a precomputed residual variance vector.

This separation is useful for marginal-likelihood or MAP tuning loops: update a
spec, build a prior, evaluate the objective, and repeat without recomputing the
residual variance scale.

## See Also

`minnesotamniwbvarm`, `minnesotainwbvarm`, `minnesotanbvarm`,
`estimateResidualVariances`, `fminsearch`, `fminunc`
