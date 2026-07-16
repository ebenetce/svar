# svar.irf

Compute impulse responses from a VAR model and impact matrix.

`svar.irf` multiplies VAR moving-average coefficient matrices by one or more
impact vectors to produce structural impulse responses.

## Syntax

```matlab
responses = svar.irf(varMdl,impact,horizon)
responses = svar.irf(varMdl,impact,horizon,Phi=Phi)
[responses,Phi] = svar.irf(___)
```

## Description

`responses = svar.irf(varMdl,impact,horizon)` computes impulse responses from
horizon 0 through `horizon` for each shock column in `impact`.

`responses = svar.irf(varMdl,impact,horizon,Phi=Phi)` uses precomputed
moving-average coefficient blocks. This is useful when evaluating many impact
matrices for the same VAR model.

`[responses,Phi] = svar.irf(___)` also returns the moving-average coefficient
blocks used in the calculation.

## Input Arguments

`varMdl` - VAR model
: Model with `NumSeries`, `P`, and `AR` properties, such as a `varm` object.

`impact` - Impact matrix
: Numeric matrix with one row per response series and one column per structural
  shock.

`horizon` - Maximum response horizon
: Nonnegative integer.

## Name-Value Arguments

`Phi` - Moving-average coefficient blocks
: Numeric array of size `NumSeries`-by-`NumSeries`-by-(`horizon` + 1). If
  omitted, `svar.irf` computes it using `svar.companionPower`.

## Output Arguments

`responses` - Impulse responses
: Numeric array of size (`horizon` + 1)-by-`NumSeries`-by-`NumShocks`.

`Phi` - Moving-average coefficient blocks
: Numeric array used to compute the responses.

## Examples

### Compute Responses for One Impact Vector

```matlab
horizon = 20;
impact = chol(EstMdl.Covariance,"lower");
responses = svar.irf(EstMdl,impact(:,1),horizon);
```

### Reuse Companion Powers

```matlab
Phi = svar.companionPower(EstMdl,horizon);
responses1 = svar.irf(EstMdl,impact1,horizon,Phi=Phi);
responses2 = svar.irf(EstMdl,impact2,horizon,Phi=Phi);
```

## See Also

`svar.companionPower`, `svar.companionMatrix`, `varm`, `varmFromCoefficients`
