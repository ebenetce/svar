# svar.companionPower

Compute moving-average blocks from powers of a VAR companion matrix.

`svar.companionPower` returns the selected top-left blocks of companion-matrix
powers used in impulse-response and forecast-error calculations.

## Syntax

```matlab
Phi = svar.companionPower(Mdl,maxLag)
```

## Description

`Phi = svar.companionPower(Mdl,maxLag)` computes the
`NumSeries`-by-`NumSeries` moving-average coefficient blocks for horizons 0
through `maxLag`. `Phi(:,:,h+1)` equals the selected block of the companion
matrix raised to power `h`.

## Input Arguments

`Mdl` - VAR model
: Model with `NumSeries`, `P`, and `AR` properties, such as a `varm` object.

`maxLag` - Maximum lag or horizon
: Nonnegative integer.

## Output Arguments

`Phi` - Moving-average coefficient blocks
: Numeric array of size `NumSeries`-by-`NumSeries`-by-(`maxLag` + 1).

## Examples

### Precompute Blocks for Impulse Responses

```matlab
horizon = 20;
Phi = svar.companionPower(EstMdl,horizon);
responses = svar.irf(EstMdl,impact,horizon,Phi=Phi);
```

## More About

### Reuse Across Shocks

The companion-power blocks depend only on the VAR model and horizon. Compute
them once per model and reuse them when evaluating multiple impact matrices or
multiple shocks.

## See Also

`svar.companionMatrix`, `svar.irf`, `varm`
