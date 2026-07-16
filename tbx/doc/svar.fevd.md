# svar.fevd

Compute forecast error variance decompositions from structural responses.

`svar.fevd` computes the share of each variable's forecast error variance
attributable to each structural shock column in an impact matrix.

## Syntax

```matlab
decomposition = svar.fevd(varMdl,impact,horizon)
decomposition = svar.fevd(varMdl,impact,horizon,Phi=Phi)
[decomposition,responses,Phi] = svar.fevd(___)
```

## Description

`decomposition = svar.fevd(varMdl,impact,horizon)` computes FEVDs from horizon
0 through `horizon`. The function first computes structural impulse responses
using `svar.irf`, squares and cumulatively sums those responses over horizons,
and normalizes by each variable's total forecast error variance.

`decomposition = svar.fevd(varMdl,impact,horizon,Phi=Phi)` uses precomputed
moving-average coefficient blocks. This is useful when evaluating multiple
impact matrices for the same VAR model.

`[decomposition,responses,Phi] = svar.fevd(___)` also returns the impulse
responses and moving-average coefficient blocks used in the calculation.

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
  omitted, `svar.fevd` computes it through `svar.irf`.

## Output Arguments

`decomposition` - Forecast error variance decomposition
: Numeric array of size (`horizon` + 1)-by-`NumSeries`-by-`NumShocks`.

`responses` - Impulse responses
: Numeric array used to compute the FEVD.

`Phi` - Moving-average coefficient blocks
: Numeric array used to compute the responses.

## Examples

### Compute FEVDs for Orthogonalized Shocks

```matlab
horizon = 20;
impact = chol(EstMdl.Covariance,"lower");
decomposition = svar.fevd(EstMdl,impact,horizon);
```

### Reuse Companion Powers

```matlab
Phi = svar.companionPower(EstMdl,horizon);
decomposition1 = svar.fevd(EstMdl,impact1,horizon,Phi=Phi);
decomposition2 = svar.fevd(EstMdl,impact2,horizon,Phi=Phi);
```

## See Also

`svar.irf`, `svar.companionPower`, `svar.companionMatrix`, `fevd`
