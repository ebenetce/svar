# marginalLikelihood

Marginal likelihood for supported Bayesian VAR priors.

## Syntax

```matlab
ml = marginalLikelihood(Mdl,Y)
ml = marginalLikelihood(Mdl,Y,X=X,Y0=Y0)
```

## Description

`ml = marginalLikelihood(Mdl,Y)` returns `p(Y | Mdl)` by exponentiating
`logMarginalLikelihood(Mdl,Y)`.

Use `logMarginalLikelihood` for numerical work and optimization. This function
is a convenience for cases where the original likelihood scale is needed.

## See Also

`logMarginalLikelihood`
