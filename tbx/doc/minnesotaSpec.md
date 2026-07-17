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
: Return optimizer starts, lower bounds, upper bounds, and field names for
  free hyperparameters. Numeric bounds use a geometric-mean start; hyperprior
  fields use their `Bounds` and `X0` properties.

`unpack`
: Write an optimizer vector back to selected hyperparameter fields, returning a
  resolved numeric copy of the spec.

`logHyperprior`
: Evaluate log-density terms from embedded `hyperprior` fields against a
  resolved candidate spec.

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

### Tune with Embedded Hyperpriors

```matlab
spec = minnesotaSpec("mniw", ...
    lambda1=hyperprior("Gamma",0.2,0.4,Bounds=[1e-4 5]), ...
    lambda3=1, ...
    lambda4=hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
    lambda5=hyperprior("Gamma",1,1,Bounds=[1e-4 50]));

psi0 = estimateResidualVariances(Y,1,Method="conditional");
Psi = arrayfun(@(x) hyperprior("InverseGamma",0.02^2,0.02^2, ...
    X0=x, Bounds=[1/100 100]*x), psi0);

PriorMdl = glp(size(Y,2),numLags,Y,Psi,Spec=spec);
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
