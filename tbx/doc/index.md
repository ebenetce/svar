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
* [minnesotaSpec](minnesotaSpec.md) - Public entry point for constructing and
  tuning Minnesota priors. Use `build` to create the concrete model object.
* [glp](glp.md) - Tune conjugate Minnesota hyperparameters by marginal
  likelihood.
* [hyperprior](hyperprior.md) - Scalar hyperprior for Minnesota-prior tuning.

Minnesota implementation reference: `minnesotaSpec` returns
[minnesotamniwSpec](minnesotamniwSpec.md),
[minnesotainwSpec](minnesotainwSpec.md), or
[minnesotanSpec](minnesotanSpec.md). Their `build` methods return
[minnesotamniwbvarm](minnesotamniwbvarm.md),
[minnesotainwbvarm](minnesotainwbvarm.md), or
[minnesotanbvarm](minnesotanbvarm.md). The shared base classes
[minnesotaBaseSpec](minnesotaBaseSpec.md) and
[minnesotabvarm](minnesotabvarm.md) are documented for maintenance reference.

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
