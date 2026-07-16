# svar.companionMatrix

Construct the companion matrix for a VAR model.

`svar.companionMatrix` stacks a VAR(`p`) model into first-order companion form.

## Syntax

```matlab
companion = svar.companionMatrix(Mdl)
```

## Description

`companion = svar.companionMatrix(Mdl)` returns the block companion matrix for
the autoregressive coefficients in `Mdl`. The first block row contains the VAR
lag coefficient matrices and the lower block rows shift lagged states forward.

## Input Arguments

`Mdl` - VAR model
: Model with `NumSeries`, `P`, and `AR` properties, such as a `varm` object.

## Output Arguments

`companion` - Companion matrix
: Numeric matrix of size (`NumSeries`*`P`)-by-(`NumSeries`*`P`).

## Examples

### Convert a VAR Model to Companion Form

```matlab
companion = svar.companionMatrix(EstMdl);
eigenvalues = eig(companion);
```

## See Also

`svar.companionPower`, `svar.irf`, `varm`
