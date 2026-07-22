# svar.estimateResidualVariances

Estimate per-series residual variances from univariate autoregressions.

`svar.estimateResidualVariances` fits an independent AR model to each response
series and returns the innovation variance estimates used as the residual
variance scale for Minnesota-prior workflows.

## Syntax

```matlab
psi = svar.estimateResidualVariances(Y,numLags)
psi = svar.estimateResidualVariances(Y,numLags,Method=method)
```

## Description

`psi = svar.estimateResidualVariances(Y,numLags)` fits an AR(`numLags`) model
to each column of `Y` using exact Gaussian maximum likelihood and returns a
1-by-`NumSeries` row vector of residual variance estimates.

`psi = svar.estimateResidualVariances(Y,numLags,Method=method)` chooses the
estimation method. Use `Method="conditional"` for conditional maximum
likelihood based on the observed lagged sample.

## Input Arguments

`Y` - Response data
: Numeric matrix, table, or timetable. Variables or columns are treated as
  separate response series.

`numLags` - Autoregressive lag order
: Positive integer.

## Name-Value Arguments

`Method` - Residual variance estimator
: `"exact"` (default) | `"conditional"`.

  `"exact"` fits univariate `arima(numLags,0,0)` models and uses the estimated
  innovation variances. `"conditional"` fits one-variable `varm` models using
  the conditional likelihood.

## Output Arguments

`psi` - Residual variance estimates
: 1-by-`NumSeries` row vector.

## Examples

### Estimate Minnesota Prior Scale

```matlab
numLags = 4;
psi = svar.estimateResidualVariances(Y,numLags);
PriorMdl = minnesotabvarm(size(Y,2),numLags,Y,Psi=psi);
```

### Use the Conditional Estimator

```matlab
psi = svar.estimateResidualVariances(Y,4,Method="conditional");
```

Use the conditional method when you want a direct noniterative estimate and do
not need exact reproducibility with workflows that use exact Gaussian
likelihood.

## More About

### Reuse in Tuning Loops

The residual variance vector is prior scale information, not a hyperparameter
candidate. Compute it once for a given data set and lag order, then reuse it
while varying `minnesotabvarm` or `svar.minnesotamniwbvarm` hyperparameters.

## See Also

`minnesotabvarm`, `svar.minnesotamniwbvarm`, `arima`, `varm`, `estimate`
