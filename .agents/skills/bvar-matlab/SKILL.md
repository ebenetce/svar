---
name: bvar-matlab
description: >-
  Build, estimate, simulate, identify, and apply Bayesian and Structural Vector
  Autoregression (BVAR / SVAR) models in MATLAB using the existing Econometrics
  Toolbox and this repo's svar functions. Use this skill
  WHENEVER the user mentions a Bayesian VAR, BVAR, SVAR, Structural VAR, a
  Minnesota/Litterman prior, a Jeffreys/diffuse or Normal-Inverse-Wishart or
  (semi)conjugate VAR prior, BVAR forecasting, or BVAR/SVAR impulse responses /
  FEVD / historical decomposition in MATLAB — even if they don't name the exact
  function. Also trigger when a user has macro/finance time series and asks for
  a shrinkage VAR, multivariate Bayesian forecasting, or structural shock
  analysis in MATLAB. Produces a runnable .m script covering the pipeline:
  prior → estimate → simulate → identify → apply.Produces a runnable .m script covering the needed
  stages of the pipeline.
---

# Bayesian / Structural VAR (BVAR / SVAR) in MATLAB

BVAR/SVAR modeling is the same handful of stages every time. Whatever the user
is after - forecast a system, run a structural shock analysis, reproduce a
study, teach the method - you build it faithfully and in a few lines by
picking the prior that matches the request and running the pipeline:

1. **Choose a prior** - the one that matches the request.
2. **Estimate** - condition it on the data, getting a posterior object.
3. **Simulate** - draw reduced-form parameters from the posterior.
4. **Identify** - turn a covariance draw into a structural `impact` matrix.
5. **Apply** - IRF, FEVD, historical decomposition, forecasts.

Not every request needs all five. A forecasting BVAR runs stages 1-3 and then
forecasts (no `impact` matrix, so stage 4 is skipped); structural work
(IRF/FEVD/HD) is what needs stage 4. Pick the stages the request actually calls
for.

The deliverable is a clean, runnable, commented `.m` file sectioned with `%%`
cells - not a chat explanation.

The unifying fact: **every prior object exposes the same `estimate` /
`simulate` / `forecast` surface.** Once stage 1 picks the object, stages 2-3
are written the same way regardless of which prior it is; stages 4-5 operate
on a `varm` plus a structural `impact` matrix. That uniformity is what keeps
the code minimal.

## Core principle: minimal code

Each stage is essentially one call. Reproduce, don't reimplement.

- **Match the model, then let the infrastructure do the rest.** Pick the prior
  that matches the request; then `estimate` / `simulate` /
  `svar.varmFromCoefficients` / `svar.irf` / `svar.fevd` /
  `historicalDecomposition` carry the analysis. Do **not** hand-roll a
  draw-to-`varm` converter, re-derive `Mu`/`V`, or write a manual IRF loop -
  those already exist.
- **Prefer natural output shapes; minimise `reshape` / `permute` / `squeeze`.**
  `svar.irf` and `svar.fevd` already return `(horizon+1) x series x shock`, so
  quantile or plot their output directly.
- **If the infrastructure is genuinely missing something**, a local function is
  fine for the analysis at hand - and flag in one line whether it is general
  enough to belong in the toolbox, so the infrastructure grows in the right
  direction. Building new prior classes, identification schemes, or apply
  functions is toolbox development, governed by the repo's `AGENTS.md`.

---

## Stage 1 - Choose the prior that matches the request

The prior encodes the modeling assumptions; pick the object whose assumptions
match what's being asked, and everything downstream is identical.

| If the analysis calls for... | Use | `estimate` returns |
|---|---|---|
| a **diffuse / Jeffreys / flat** prior (OLS-equivalent recursive SVARs; e.g. Kilian 2009, textbook VARs) | `diffusebvarm(m,p)` | `conjugatebvarm` (analytic) |
| a **Minnesota / Litterman** shrinkage prior (forecasting, large systems; e.g. GLP 2012, Banbura-Giannone-Reichlin 2010) | `minnesotabvarm(m,p,Y,...)` (+ `glp` to tune) | analytic / Gibbs by `Type` |
| **Uhlig's weak** improper NIW | `weakbvarm(m,p)` | `conjugatebvarm` (analytic) |
| a **uniform prior over impact rotations** (set/sign-ID designs) | `uniformirbvarm(m,p)` | `conjugatebvarm` (analytic) |
| an **independent Normal-Wishart** (asymmetric own/cross shrinkage) | `minnesotabvarm(...,Type="inw")` or `semiconjugatebvarm(m,p,...)` | `empiricalbvarm` (Gibbs) |
| **normal coefficients, fixed Sigma** | `minnesotabvarm(...,Type="normal")` or `normalbvarm(m,p,...)` | `normalbvarm` (analytic) |
| a **fully specified Matrix-Normal-IW** you set yourself | `conjugatebvarm(m,p,...)` with `Mu,V,Omega,DoF` | `conjugatebvarm` (analytic) |

`m = numseries`, `p = numlags`. All accept `SeriesNames`, `IncludeConstant`,
`IncludeTrend`, `NumPredictors`. (`bayesvarm` is only a Minnesota shortcut - use
`minnesotabvarm` unless the user asks for it.)

**No prior at all - a classical VAR?** If the request is a plain OLS VAR (point
estimates, bootstrap bands, no Bayesian prior), skip the BVAR object and fit
`varm` / `estimate` directly. Only stages 1-3 change: stages 4-5 are
prior-agnostic - `svar.irf` / `svar.fevd` /
`historicalDecomposition` take **any** `varm`. Reach for `diffusebvarm` instead
when you want the Bayesian counterpart: the same OLS point estimates plus a
posterior and analytic credible bands.

Two families need a word more.

**Reduced-form improper priors** (`weakbvarm`, `uniformirbvarm`) are data-free -
pass `Y` at `estimate`:
```matlab
prior = weakbvarm(m, p);                              % Uhlig (2005) weak improper NIW
prior = uniformirbvarm(m, p, DeterminantShift=-1);    % uniform over impact rotations
```

**The Minnesota family** - the `minnesotabvarm` front door takes `Y` up front
and returns a complete prior (residual-variance scale and any dummy
observations already resolved, sample stored on the object):
```matlab
prior = minnesotabvarm(m, p, Y);                        % conjugate MNIW (default)
prior = minnesotabvarm(m, p, Y, lambda4=1, lambda5=1);  % + sum-of-coeff & dummy-initial-obs
prior = minnesotabvarm(m, p, Y, Type="inw", lambda2=0.5);
```
`Type` is `"mniw"` (default), `"inw"`, or `"normal"` (aliases ignore
hyphens/underscores/spaces). Hyperparameters: `lambda1` overall tightness,
`lambda3` lag decay, `lambda4`/`lambda5` sum-of-coefficients &
dummy-initial-observation (`mniw`), `lambda2` cross-tightness (`inw`/`normal`),
`Vc`, `PriorMean` (`0` for differenced series), `Psi` residual-variance scale.
The family's own constructor validates options - passing `lambda2` to `"mniw"`
errors. Full detail: `references/priors.md`.

**Tune the Minnesota tightnesses by marginal likelihood** (conjugate MNIW) with
`glp` (Giannone-Lenza-Primiceri) - the way to reproduce GLP-style hierarchical
BVARs. A hyperparameter is *free* if you pass a `hyperprior` or `[lo hi]`
bounds, *fixed* if you pass a scalar:
```matlab
psi = svar.estimateResidualVariances(Y, p, Method="conditional");   % 1-by-m scale
[prior, info] = glp(m, p, Y, ...
    lambda1 = hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5]), ...   % free, Gamma prior
    lambda4 = [1e-4 50], ...                                        % free, flat prior
    lambda5 = 1, ...                                                % fixed
    Psi     = psi);
```
`glp` returns the prior at the marginal-likelihood maximiser. To integrate
*over* the hyperparameters (the full hierarchical treatment), ask for a third
output and a Metropolis chain:
```matlab
[prior, info, chain] = glp(m, p, Y, lambda1=..., Psi=psi, ...
    NumDraws=10000, BurnIn=10000, ProposalScale=0.8);
% chain.Coefficients / chain.Sigma already integrate over the hyperparameters
```
Feed `chain.Coefficients(:,:,d)` / `chain.Sigma(:,:,d)` straight into stage 5;
they skip stages 2-3. See `references/pipeline.md`.

**If no existing prior matches the request**, a bespoke prior class is toolbox
development, not analysis-script code - see `AGENTS.md` (subclassing, hand-built
`Mu`/`V`, Sims-Zha dummies). Don't inline a new prior class into the script.

---

## Stage 2 - Estimate

Condition the prior on the data. Identical call for every prior object:

```matlab
[Posterior, Summary] = estimate(prior, Y, Display="off");
```

**Object-state rule.** Priors that already carry the sample -
`minnesotabvarm` and the `glp` output - are estimated with **no** `Y` argument
(`estimate(prior)`); passing a *different* `Y` errors, because their scale and
dummies derive from the stored sample. Data-free priors - everything except
`minnesotabvarm` and the `glp` output - **require** `Y`. The posterior class
depends on the family:

| Prior | Posterior class | Draws come from |
|---|---|---|
| `conjugatebvarm`, `diffusebvarm`, `weakbvarm`, `uniformirbvarm`, Minnesota `mniw` | `conjugatebvarm` | analytic NIW |
| `normalbvarm`, Minnesota `normal` | `normalbvarm` | analytic (Sigma fixed) |
| `semiconjugatebvarm`, Minnesota `inw` | `empiricalbvarm` | Gibbs sampler |

Set `rng(...)` before estimating any **sampler-based** (`inw` / semiconjugate)
model. Optional `X=` (exogenous predictors) and `Y0=` (presample) are accepted.

---

## Stage 3 - Simulate

Draw reduced-form parameters for structural work. Same call shape everywhere:

```matlab
rng(1);
[CoeffDraws, SigmaDraws] = simulate(Posterior, NumDraws=1000);
% CoeffDraws: numcoeff-by-NumDraws ;  SigmaDraws: m-by-m-by-NumDraws
```

The same **object-state rule** applies: call `simulate(Posterior, NumDraws=n)`
when the object embeds the sample (a posterior from `estimate`, a
`minnesotabvarm`, the `glp` output); call `simulate(prior, Y, NumDraws=n)` only
when the object is a data-free prior that still needs conditioning. Never pass
`Y` to something that already carries it - see `references/priors.md` section
6.

For semiconjugate / `inw` models `estimate` returns an `empiricalbvarm`
summary; when you need raw draws, prefer `simulate(prior, Y, NumDraws=n)`.

---

## Stage 4 - Identify

Identification turns each reduced-form covariance draw `Sigma` into a
structural **impact matrix** `impact` (`numSeries x numShocks`) - the
contemporaneous map from structural shocks to responses.

The toolbox currently ships **recursive (Cholesky)** identification, which is
one line:

```matlab
impact = chol(Sigma, "lower");    % recursive ID; the SeriesNames order IS the assumption
```

**Always comment the variable ordering** wherever Cholesky is used - it *is*
the identifying assumption. To flip a shock's sign, negate that shock's
`impact` column (VARs are linear).

> Sign- and zero-restriction identification (Haar rotations,
> Arias-Rubio-Ramirez-Waggoner) is on the roadmap and will land as a dedicated
> identification API. Until then, supply `impact` yourself as above; a
> sign/zero routine you write for one analysis is local glue, and adding a
> scheme to the toolbox is an `AGENTS.md` task.

---

## Stage 5 - Apply

All apply functions take a `varm` object plus the `impact` matrix. Bridge from
a BVAR draw or posterior to a `varm` with `svar.varmFromCoefficients` /
`bvar2var` - never assemble a `varm` by hand:

```matlab
Mdl = svar.varmFromCoefficients(Posterior, reshape(CoeffDraws(:,d), [], m), SigmaDraws(:,:,d)); % per draw
Mdl = bvar2var(Posterior);                                                      % whole model
```

`simulate` returns each coefficient draw as a `numcoeff x 1` vector, so reshape
it to `m*k x m` (`reshape(CoeffDraws(:,d), [], m)`) before
`svar.varmFromCoefficients`. A `glp` chain already stores
`chain.Coefficients(:,:,d)` in matrix form - pass it straight through.

`svar.companionPower` is the expensive, impact-independent piece - compute it
**once per `varm`** and reuse it across `svar.irf` / `svar.fevd`:

```matlab
H   = 20;
Phi = svar.companionPower(Mdl, H);              % reduced-form MA blocks, reused below
ir  = svar.irf(Mdl,  impact, H, Phi=Phi);       % (H+1) x series x shock
fe  = svar.fevd(Mdl, impact, H, Phi=Phi);       % (H+1) x series x shock shares
```

`svar.irf` / `svar.fevd` already return `horizon x response x shock`, so
quantile or plot them directly - **no `permute`/`reshape`**.

`historicalDecomposition` takes a **BVAR object or a `varm`**, a **square**
`impact`, and `Y`, and returns a struct (`Contributions` `T-p x series x
shock`, `Total`, `StructuralShocks`, ...). Passing the BVAR posterior directly
is simplest:

```matlab
impact = chol(bvar2var(Posterior).Covariance, "lower");   % comment the ordering
HD     = historicalDecomposition(Posterior, impact, Y);   % or a varm from bvar2var
```

The canonical credible-band loop (simulate -> per-draw
`svar.varmFromCoefficients` -> identify -> `svar.irf`/`svar.fevd` -> quantiles)
and the FEVD/HD plotting patterns are in `references/pipeline.md`.

**Forecasts** use `forecast` / `simsmooth` directly on the posterior object:
`forecast(Posterior, fh, Y)` for mean/variance bands, or `simsmooth` with
NaN-padded data for full predictive paths (see `references/pipeline.md`).
Conditional forecasts are on the roadmap; until then, condition via
`simsmooth` with the known future values pinned in the NaN-padded sample.

---

## SSVS - what's available

SSVS is **not** available for the BVAR (`*bvarm`) objects.
It *is* available for Bayesian linear regression via `mixsemiconjugateblm` (and
custom priors via `customblm`). For a VAR: either approximate with heavy
shrinkage in a semiconjugate / `inw` Minnesota prior, or estimate
equation-by-equation with `mixsemiconjugateblm` and assemble the system. State
which you're doing; don't imply the BVAR objects support SSVS.

## Reference files
- `references/priors.md` - Stage 1 detail: the prior constructors, the
  `minnesotabvarm` / `glp` / `hyperprior` tuning API, the hyperparameter table,
  and the object-state rule. Read before choosing or tuning a prior.
- `references/pipeline.md` - Stages 2-5 detail: lag selection, forecast
  intervals / fan charts, the per-draw and hierarchical IRF/FEVD credible-band
  loops, FEVD plotting, and historical decomposition.

## Output conventions
- Keep it minimal: build the requested model, then call the existing
  functions. No reimplementation of what ships. Any local helper comes with a
  one-line "should this live in `tbx/svar`?" note; building new toolbox classes
  is an `AGENTS.md` task.
- Choose the prior by fit to the request.
- Set `rng` before any sampler-based (`inw` / semiconjugate) estimation or
  simulation, and before any `glp` Metropolis chain.
- Respect the object-state rule: pass `Y` only to data-free priors.
- Comment the variable ordering wherever Cholesky identification is used.
- Prefer natural output shapes; avoid `reshape`/`permute`/`squeeze` unless a
  function's documented output genuinely requires it.
- Section scripts with `%%` cell dividers so they run cell-by-cell.
