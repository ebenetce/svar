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
spec = minnesotaSpec(Method,Name=Value)
```

### Description

`spec = minnesotaSpec` creates a `minnesotamniwSpec` with default Minnesota
hyperparameters.

`spec = minnesotaSpec(Method)` chooses `"mniw"`, `"inw"`, or `"normal"`.
Hyphens, underscores, and spaces are ignored when matching method aliases.

`spec = minnesotaSpec(Name=Value)` sets default MNIW hyperparameter fields.

`spec = minnesotaSpec(Method,Name=Value)` chooses a specification family and
sets one or more hyperparameter fields.

## Input Arguments

`Method` - Specification family
: `"mniw"` (default) | `"conjugate"` | `"matrixnormal"` | `"inw"` |
  `"independent"` | `"kadiyala"` | `"kadiyalakarlsson"` |
  `"semiconjugate"` | `"normal"` | `"litterman"` | `"fixed"` |
  `"fixedsigma"`.

## Name-Value Arguments

The available name-value arguments depend on `Method`.

`lambda1` - Overall tightness
: `0.2` (default) | positive scalar | `[lower upper]` bounds |
  `hyperprior`.

`lambda2` - Cross-variable relative tightness
: Available for `"inw"` and `"normal"` specs. `0.5` (default) | positive
  scalar | `[lower upper]` bounds | `hyperprior`.

`lambda3` - Lag-decay exponent
: `1` (default) | nonnegative scalar | `[lower upper]` bounds |
  `hyperprior`.

`lambda4` - Sum-of-coefficients dummy tightness
: Available for `"mniw"` specs. `Inf` (default) | positive scalar |
  `[lower upper]` bounds | `hyperprior`.

`lambda5` - Dummy-initial-observation tightness
: Available for `"mniw"` specs. `Inf` (default) | positive scalar |
  `[lower upper]` bounds | `hyperprior`.

`Vc` - Prior variance for deterministic terms
: `1e4` (default) | positive scalar.

`PriorMean` - Own first-lag prior mean
: `[]` (default) | numeric row vector. An empty value means the built model
  uses a row vector of ones.

## Output Arguments

`spec` - Minnesota specification
: `minnesotamniwSpec`, `minnesotainwSpec`, or `minnesotanSpec`.

## Object Functions

`build`
: Materialize a concrete Minnesota model from the spec and a residual variance
  vector.

`pack`
: Return optimizer starts, lower bounds, upper bounds, and field names for free
  hyperparameters.

`unpack`
: Write optimizer values back to selected hyperparameter fields.

`freeFields`
: Return the names of free hyperparameters.

`logHyperprior`
: Available on MNIW specs. Evaluate log-density terms from embedded
  `hyperprior` fields.

## Examples

### Build a Prior from a Spec

```matlab
numLags = 4;
psi = estimateResidualVariances(Y,numLags,Method="conditional");

spec = minnesotaSpec(lambda1=0.25,lambda4=10,lambda5=5);
PriorMdl = spec.build(size(Y,2),numLags,psi, ...
    SeriesNames=["Output" "Prices" "Rate"]);
```

### Select a Nonconjugate Specification

```matlab
spec = minnesotaSpec("inw",lambda1=0.2,lambda2=0.4,lambda3=1);
PriorMdl = spec.build(size(Y,2),4,psi);
```

### Pack and Unpack Free Parameters

```matlab
spec = minnesotaSpec("mniw",lambda1=[0.01 1],lambda4=[1 50]);
[theta0,lb,ub,names] = spec.pack();

candidate = spec.unpack(theta0,names);
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

### Specification Families

The MNIW family builds `minnesotamniwbvarm` objects and supports an analytic
marginal likelihood. The INW family builds `minnesotainwbvarm` objects and
supports free cross-variable shrinkage, $$\lambda_2$$, with a simulation-based
posterior. The Normal family builds `minnesotanbvarm` objects with fixed
innovations covariance.

### Scalar-or-Bounds Convention

For tunable hyperparameter fields, a scalar fixes the value, a two-element
vector `[lower upper]` makes the field free with flat bounds, and a `hyperprior`
object makes the field free with an added log-density term.

## See Also

`minnesotamniwSpec`, `minnesotainwSpec`, `minnesotanSpec`,
`minnesotamniwbvarm`, `minnesotainwbvarm`, `minnesotanbvarm`, `hyperprior`,
`glp`

