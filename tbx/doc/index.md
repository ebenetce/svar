# SVAR Toolbox

Structural VAR workflows for MATLAB built on top of Econometrics Toolbox VAR
and Bayesian VAR model objects.

This toolbox adds reduced-form Bayesian VAR priors, Minnesota-prior
hyperparameter utilities, and structural impulse-response helpers for research
workflows that start from `conjugatebvarm` or `varm` models.

## Bayesian VAR Priors

* [weakbvarm](weakbvarm.md) - Weak improper Bayesian VAR prior for
  Uhlig-style reduced-form workflows.
* [uniformirbvarm](uniformirbvarm.md) - Reduced-form prior induced by uniform
  impulse-response parameters.
* [minnesotabvarm](minnesotabvarm.md) - Create a Minnesota prior by selecting
  a supported prior family with the `Type` name-value argument.
* [svar.minnesotamniwbvarm](minnesotamniwbvarm.md) - Conjugate
  Matrix-Normal-Inverse-Wishart Minnesota prior.
* [svar.minnesotainwbvarm](minnesotainwbvarm.md) - Semiconjugate independent
  Normal-Wishart Minnesota prior.
* [svar.minnesotanbvarm](minnesotanbvarm.md) - Normal Minnesota prior with
  fixed innovations covariance.
* [glp](glp.md) - Tune conjugate Minnesota hyperparameters by marginal
  likelihood.
* [logMarginalLikelihood](logMarginalLikelihood.md) - Compute analytic log
  evidence for supported Bayesian VAR priors.
* [marginalLikelihood](marginalLikelihood.md) - Return marginal likelihood on
  the original scale.
* [hyperprior](hyperprior.md) - Scalar hyperprior for Minnesota-prior tuning.

`minnesotabvarm` is the recommended front door. It dispatches to
`svar.minnesotamniwbvarm`, `svar.minnesotainwbvarm`, or
`svar.minnesotanbvarm`.

## Reduced-Form Utilities

* [svar.estimateResidualVariances](estimateResidualVariances.md) - Estimate
  the residual variance scale used by Minnesota priors.
* [svar.varmFromCoefficients](varmFromCoefficients.md) - Reconstruct a `varm`
  model from a Bayesian VAR coefficient draw.
* [bvar2var](bvar2var.md) - Convert a Bayesian VAR model to a `varm` model.
* [svar.logDetPD](logDetPD.md) - Compute the log determinant of a
  positive-definite matrix using a Cholesky factorization.

## Structural VAR Package Functions

* [svar.irf](svar.irf.md) - Compute structural impulse responses from a VAR
  model and impact matrix.
* [svar.fevd](svar.fevd.md) - Compute forecast error variance decompositions
  from structural impulse responses.
* [historicalDecomposition](historicalDecomposition.md) - Decompose a realized
  history into structural shock contributions.
* [svar.companionMatrix](svar.companionMatrix.md) - Build the VAR companion
  matrix.
* [svar.companionPower](svar.companionPower.md) - Compute moving-average blocks
  from companion-matrix powers.

---
