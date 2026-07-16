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
* [minnesotabvarm](minnesotabvarm.md) - Conjugate Minnesota prior with optional
  dummy-observation priors.
* [minnesotaSpec](minnesotaSpec.md) - Hyperparameter recipe for constructing
  and tuning Minnesota priors.

## Reduced-Form Utilities

* [estimateResidualVariances](estimateResidualVariances.md) - Estimate the
  residual variance scale used by Minnesota priors.
* [varmFromCoefficients](varmFromCoefficients.md) - Reconstruct a `varm` model
  from a Bayesian VAR coefficient draw.
* [bvar2var](bvar2var.md) - Convert a Bayesian VAR model to a `varm` model.
* [logDetPD](logDetPD.md) - Compute the log determinant of a positive-definite
  matrix using a Cholesky factorization.

## Structural VAR Package Functions

* [svar.irf](svar.irf.md) - Compute structural impulse responses from a VAR
  model and impact matrix.
* [svar.fevd](svar.fevd.md) - Compute forecast error variance decompositions
  from structural impulse responses.
* [svar.companionMatrix](svar.companionMatrix.md) - Build the VAR companion
  matrix.
* [svar.companionPower](svar.companionPower.md) - Compute moving-average blocks
  from companion-matrix powers.

---
