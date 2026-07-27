# svar.seasonalDummies

Create seasonal indicator variables.

`svar.seasonalDummies` creates a numeric seasonal-dummy matrix for use as
exogenous predictors in VAR, BVAR, and SVAR workflows.

## Syntax

```matlab
X = svar.seasonalDummies(numObs,period)
X = svar.seasonalDummies(...,DropLast=tf)
X = svar.seasonalDummies(...,Reference=ref)
```

## Description

`X = svar.seasonalDummies(numObs,period)` creates a matrix with `numObs`
rows and repeating seasonal categories `1:period`. The first observation is
assigned to category 1. By default, the last category is omitted so `X` can be
used with a model intercept.

`X = svar.seasonalDummies(...,DropLast=tf)` controls whether one category is
omitted. Set `DropLast=false` to return all `period` columns.

`X = svar.seasonalDummies(...,Reference=ref)` specifies the omitted category
when `DropLast=true`. The default is `period`.

## Input Arguments

`numObs` - Number of observations
: Positive integer.

`period` - Seasonal period
: Positive integer. Use `12` for monthly data or `4` for quarterly data.

## Name-Value Arguments

`DropLast` - Drop reference category
: Logical scalar. The default is `true`.

`Reference` - Reference category
: Positive integer less than or equal to `period`. The default is `period`.

## Output Arguments

`X` - Seasonal indicators
: Numeric matrix. The size is `numObs`-by-(`period` - 1) when `DropLast=true`
and `numObs`-by-`period` otherwise.

## Examples

### Monthly Dummies With an Intercept

Create monthly dummy variables for a sample that starts in February. Since the
first observation is category 1, the default omitted category 12 corresponds to
the following January.

```matlab
X = svar.seasonalDummies(36,12);
```

### Quarterly Dummies Without an Intercept

```matlab
X = svar.seasonalDummies(20,4,DropLast=false);
```

## See Also

`diffusebvarm`, `estimate`, `varm`
