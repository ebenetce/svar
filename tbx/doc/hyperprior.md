# hyperprior

Scalar hyperprior for Minnesota-prior hyperparameter tuning.

The `hyperprior` class stores one probability distribution for one scalar
hyperparameter. It provides optimizer bounds, an initial value, and a stable
log-density used by `glp`.

## Creation

### Syntax

```matlab
hp = hyperprior(distribution,mode,sd)
hp = hyperprior(distribution,mode,sd,Bounds=bounds)
hp = hyperprior(distribution,mode,sd,X0=x0)
```

### Description

`hp = hyperprior(distribution,mode,sd)` creates a hyperprior from an
interpretable mode and standard deviation. The object converts those moments to
the native distribution parameters once at construction time.

`hp = hyperprior(distribution,mode,sd,Bounds=bounds)` sets the optimizer box.
If `Bounds` is omitted, quantile bounds are computed from the distribution.

`hp = hyperprior(distribution,mode,sd,X0=x0)` sets the optimizer starting
value. If `X0` is omitted, the distribution mode is used when valid; otherwise
the median is used.

## Input Arguments

`distribution` - Distribution family
: `"Gamma"` | `"InverseGamma"` | `"Beta"`.

`mode` - Distribution mode
: Positive scalar.

`sd` - Distribution standard deviation
: Positive scalar.

## Name-Value Arguments

`Bounds` - Optimizer bounds
: Two-element vector `[lo hi]` inside the distribution support.

`X0` - Optimizer starting point
: Finite scalar inside `Bounds` and inside the distribution support.

## Properties

`Distribution` - Distribution family
: String scalar.

`Params` - Native distribution parameters
: Two-element numeric vector. Gamma and inverse-gamma use shape and scale;
  beta uses alpha and beta.

`Bounds` - Optimizer bounds
: Two-element numeric vector.

`X0` - Optimizer starting value
: Numeric scalar.

## Object Functions

`logpdf`
: Evaluate the log density. Values outside the support return `-Inf` so an
  optimizer can probe invalid regions without throwing an error.

`quantileBounds`
: Return a quantile-based optimizer box.

`initialValue`
: Return the default optimizer starting value.

`modeOf`
: Return the mode implied by `Params`.

`meanOf`
: Return the mean implied by `Params`, when it exists.

## Examples

### Put a Gamma Hyperprior on Overall Tightness

```matlab
hp = hyperprior("Gamma",0.2,0.4,Bounds=[1e-4 5]);
PriorMdl = glp(size(Y,2),4,Y,lambda1=hp);
```

### Use Inverse-Gamma Priors for Residual Variances

```matlab
psi0 = svar.estimateResidualVariances(Y,1,Method="conditional");
Psi = arrayfun(@(x) hyperprior("InverseGamma",0.02^2,0.02^2, ...
    X0=x,Bounds=[1/100 100]*x), psi0);

PriorMdl = glp(size(Y,2),4,Y,Psi=Psi);
```

### Evaluate a Log Density

```matlab
hp = hyperprior("Beta",0.5,0.1);
lp = hp.logpdf(0.4);
```

## More About

### Moment Conversion

The constructor accepts an interpretable pair, mode and standard deviation, and
stores native distribution parameters. For a gamma distribution with shape
$$k$$ and scale $$\theta$$,

$$
\operatorname{mode}(X) = (k - 1)\theta,\qquad
\operatorname{Var}(X) = k\theta^2.
$$

Beta and inverse-gamma conversions are solved once with a scalar root finder.
The resulting object does not repeat that solve during optimization.

### Role in Tuning

Assigning a `hyperprior` to a `glp` name-value argument makes that field free.
The field's `Bounds` and `X0` participate in the optimizer vector, and
`logpdf` contributes a log-prior term to the objective.

## See Also

`glp`, `svar.minnesotamniwbvarm`, `logMarginalLikelihood`
