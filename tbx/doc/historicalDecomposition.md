# historicalDecomposition

Decompose a realized history into structural shock contributions.

`historicalDecomposition` computes the cumulative contribution of each
structural shock to each observed response series in a reduced-form Bayesian
VAR workflow.

## Syntax

```matlab
HD = historicalDecomposition(Mdl,Impact,Y)
HD = historicalDecomposition(Mdl,Impact,Y,Name=Value)
```

## Description

`HD = historicalDecomposition(Mdl,Impact,Y)` infers reduced-form residuals from
`Y`, maps them into structural shocks with `Impact`, and returns time-varying
contributions for each variable and shock. `Mdl` must be a Bayesian VAR model
accepted by `bvar2var`.

`HD = historicalDecomposition(Mdl,Impact,Y,Name=Value)` specifies presample
data, predictor data, precomputed residuals, or subsets of variables and
shocks to report.

## Input Arguments

`Mdl` - Bayesian VAR model
: Model object accepted by `bvar2var`.

`Impact` - Structural impact matrix
: Numeric square matrix with one row per response series and one column per
  structural shock.

`Y` - Response data
: Numeric matrix used by `infer` to compute reduced-form residuals. You can
  omit `Y` when you provide `Residuals`.

## Name-Value Arguments

`Y0` - Presample responses
: Numeric matrix passed to `infer`.

`X` - Predictor data
: Numeric matrix passed to `infer`.

`Residuals` - Precomputed reduced-form residuals
: Numeric matrix with one column per response series. When you specify
  `Residuals`, the function does not call `infer`.

`VariableIndices` - Variables to report
: Positive integer vector. The default reports all variables.

`ShockIndices` - Shocks to report
: Positive integer vector. The default reports all shocks.

## Output Arguments

`HD` - Historical decomposition
: Structure containing the model, impact matrix, residuals, structural shocks,
  full contribution array, selected contribution array, selected totals,
  selected indices, variable names, and shock names.

## Examples

### Decompose Posterior Contributions

```matlab
numLags = 4;
PriorMdl = minnesotabvarm(size(Y,2),numLags,Y);
PosteriorMdl = estimate(PriorMdl);
varMdl = bvar2var(PosteriorMdl);
Impact = chol(varMdl.Covariance,"lower");

HD = historicalDecomposition(PosteriorMdl,Impact,Y);
```

### Report a Subset of Variables and Shocks

```matlab
HD = historicalDecomposition(PosteriorMdl,Impact,Y, ...
    VariableIndices=1:2,ShockIndices=1);
```

## More About

### Contribution Array

`HD.FullContributions` is a `T`-by-`NumSeries`-by-`NumSeries` numeric array.
The first dimension is time, the second is the response variable, and the
third is the structural shock. `HD.Contributions` contains the selected subset
requested by `VariableIndices` and `ShockIndices`, and `HD.Total` sums the
selected shocks for each selected variable.

## See Also

`bvar2var`, `svar.irf`, `infer`, `varm`
