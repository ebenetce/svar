# varmFromCoefficients

Reconstruct a `varm` model from one Bayesian VAR coefficient draw.

`varmFromCoefficients` converts a coefficient matrix and covariance matrix from
a Bayesian VAR draw into a concrete `varm` model that can be used with
Econometrics Toolbox functions such as `forecast` and `irf`.

## Syntax

```matlab
Mdl = varmFromCoefficients(template,coefficients,covariance)
```

## Description

`Mdl = varmFromCoefficients(template,coefficients,covariance)` builds a `varm`
model using `template` for the model layout and using the supplied coefficient
and covariance draw for the numeric parameters.

The function is prior-agnostic. `template` can be a `conjugatebvarm` object or a
subclass such as `minnesotamniwbvarm`, prior or posterior. The function uses
only layout properties such as `NumSeries`, `P`, `SeriesNames`,
`IncludeConstant`, `IncludeTrend`, and `NumPredictors`.

## Input Arguments

`template` - Bayesian VAR template
: `bvar` object, such as `conjugatebvarm`, `minnesotamniwbvarm`,
  `weakbvarm`, or `uniformirbvarm`.

`coefficients` - Coefficient draw
: Numeric matrix with one row per equation coefficient and one column per
  response series.

  Rows must follow the `conjugatebvarm` coefficient layout: lag coefficient
  blocks first, then the constant row if present, then the trend row if
  present. Models with exogenous predictors are not supported by this
  reconstruction function.

`covariance` - Innovations covariance draw
: Numeric square matrix.

## Output Arguments

`Mdl` - Reconstructed VAR model
: `varm` model object.

## Examples

### Convert One Posterior Draw to `varm`

```matlab
numDraws = 1;
[Coeff,Sigma] = simulate(PosteriorMdl,NumDraws=numDraws);

drawMdl = varmFromCoefficients(PosteriorMdl,Coeff(:,:,1),Sigma(:,:,1));
```

### Compute an Impulse Response from a Draw

```matlab
drawMdl = varmFromCoefficients(PosteriorMdl,Coeff(:,:,draw),Sigma(:,:,draw));
drawIRF = svar.irf(drawMdl,impact,horizon);
```

## More About

### Predictor Limitation

`varmFromCoefficients` supports lag coefficients, constants, and trends. It
throws an error when `template.NumPredictors > 0` because the reconstruction
path does not currently map exogenous-predictor coefficients into a `varm`
object.

## See Also

`varm`, `simulate`, `svar.irf`, `minnesotamniwbvarm`, `conjugatebvarm`
