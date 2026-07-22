# hyperprior

Scalar hyperprior for Minnesota-prior hyperparameter tuning.

The `hyperprior` class stores one probability distribution for one scalar
hyperparameter. It provides optimizer bounds, an initial value, and a stable
log-density used by `glp` and `minnesotaSpec` tuning workflows.

## Creation

### Syntax

```matlab
hp = hyperprior(distribution,mode,sd)
hp = hyperprior(distribution,mode,sd,Bounds=bounds)
hp = hyperprior(distribution,mode,sd,X0=x0)
hp = hyperprior(distribution,p1,p2,Parameterization="native",Bounds=bounds)
```

### Description

`hp = hyperprior(distribution,mode,sd)` creates a hyperprior from an
interpretable mode and standard deviation. The object converts those moments to
the native distribution parameters once at construction time.

`hp = hyperprior(distribution,mode,sd,Bounds=bounds)` sets the optimizer box.
If `Bounds` is omitted, quantile bounds are computed from the distribution.

`hp = hyperprior(distribution,mode,sd,X0=x0)` sets the optimizer starting value.
If `X0` is omitted, the distribution mode is used when valid; otherwise the
median is used.

`hp = hyperprior(distribution,p1,p2,Parameterization="native")` takes the two
positional inputs as the native parameters themselves, skipping the moment
conversion. Use this for densities that no `(mode, sd)` pair identifies —
the moment solve requires a finite standard deviation, so it cannot reach a
gamma with shape below 1 or an inverse-gamma with shape below 2. Pass `Bounds`
as well: the default quantile box is computed with `gaminv`/`betainv`, which
does not converge in the far tail of such a diffuse density, and an unusable
default raises `hyperprior:unusableDefaultBounds` rather than being stored.

## Input Arguments

`distribution` - Distribution family
: `"Gamma"` | `"InverseGamma"` | `"Beta"`.

`mode` - Distribution mode
: Positive scalar. Read as the first native parameter when
  `Parameterization="native"`.

`sd` - Distribution standard deviation
: Positive scalar. Read as the second native parameter when
  `Parameterization="native"`.

## Name-Value Arguments

`Parameterization` - Meaning of the positional inputs
: `"moments"` (default) | `"native"`. With `"moments"` the inputs are a mode
  and a standard deviation; with `"native"` they are `Params` directly —
  gamma and inverse-gamma shape and scale, or beta alpha and beta.

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
spec = minnesotaSpec("mniw",lambda1=hp);
```

### Use Inverse-Gamma Priors for Residual Variances

The hyperprior on the residual-variance scale in Giannone, Lenza and Primiceri
(2012) is an inverse-gamma with shape and scale both equal to `0.02^2`. It has
neither a mean nor a variance, so it must be given natively.

```matlab
psi0 = estimateResidualVariances(Y,1,Method="conditional");
Psi = arrayfun(@(x) hyperprior("InverseGamma",0.02^2,0.02^2, ...
    Parameterization="native", X0=x, Bounds=[x/100 x*100]), psi0);

PriorMdl = glp(size(Y,2),5,Y,Psi=Psi);
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

Assigning a `hyperprior` to a Minnesota spec field makes that field free. The
field's `Bounds` and `X0` participate in the optimizer vector, and `logpdf`
contributes a log-prior term to the objective used by `glp`.

## See Also

`glp`, `minnesotaSpec`, `minnesotamniwSpec`

