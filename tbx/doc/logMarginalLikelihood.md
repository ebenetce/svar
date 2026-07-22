# logMarginalLikelihood

Analytic log marginal likelihood for supported Bayesian VAR priors.

## Syntax

```matlab
logML = logMarginalLikelihood(Mdl)
logML = logMarginalLikelihood(Mdl,Y)
logML = logMarginalLikelihood(Mdl,Y,X=X,Y0=Y0)
```

## Description

`logML = logMarginalLikelihood(Mdl)` returns the analytic log evidence using
the response sample stored on `Mdl`. This syntax is available for prior objects
that carry a sample, such as `svar.minnesotamniwbvarm`.

`logML = logMarginalLikelihood(Mdl,Y)` returns `log p(Y | Mdl)` for proper
`conjugatebvarm` and fixed-Sigma `normalbvarm` prior objects.

`logML = logMarginalLikelihood(Mdl,Y,X=X,Y0=Y0)` supplies exogenous predictors
or presample responses to the sufficient-statistics calculation.

## Supported Models

Supported:

* `conjugatebvarm`
* `normalbvarm`
* `svar.minnesotamniwbvarm`
* `svar.minnesotanbvarm`

Unsupported models throw explicit errors:

* `semiconjugatebvarm` and `svar.minnesotainwbvarm` have no closed-form
  marginal likelihood in this toolbox.
* `diffusebvarm`, `weakbvarm`, and `uniformirbvarm` are improper-prior
  workflows, so their evidence is not comparable without an imposed
  convention.

## Examples

```matlab
psi = svar.estimateResidualVariances(Y,4,Method="conditional");
PriorMdl = svar.minnesotamniwbvarm(size(Y,2),4,Y,Psi=psi);

logML = logMarginalLikelihood(PriorMdl);
```

## See Also

`marginalLikelihood`, `conjugatebvarm`, `normalbvarm`,
`svar.minnesotamniwbvarm`, `svar.minnesotanbvarm`
