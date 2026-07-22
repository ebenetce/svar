# bvar2var

Convert a Bayesian VAR model to a `varm` model.

`bvar2var` maps common reduced-form VAR properties from a Bayesian VAR object
to an Econometrics Toolbox `varm` object.

## Syntax

```matlab
varMdl = bvar2var(bvarMdl)
```

## Description

`varMdl = bvar2var(bvarMdl)` converts `bvarMdl` to a `varm` model by copying
properties such as autoregressive coefficients, constant, trend, beta,
covariance, and series names.

## Input Arguments

`bvarMdl` - Bayesian VAR model
: Supported `bvar` object.

## Output Arguments

`varMdl` - VAR model
: `varm` model object.

## Examples

### Convert a Conjugate Bayesian VAR

```matlab
bvarMdl = conjugatebvarm(3,2);
varMdl = bvar2var(bvarMdl);
```

## See Also

`varm`, `conjugatebvarm`, `svar.varmFromCoefficients`
