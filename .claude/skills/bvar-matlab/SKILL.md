---
name: bvar-matlab
description: >-
  Build, estimate, simulate, identify, and apply Bayesian and Structural Vector
  Autoregression (BVAR / SVAR) models in MATLAB with the svar toolbox — an
  extension of the Econometrics Toolbox. Every BVAR follows the same five-stage
  pipeline: choose a prior, estimate, simulate, identify, apply. Covers the
  built-in prior objects (conjugatebvarm, semiconjugatebvarm, normalbvarm,
  diffusebvarm) AND this toolbox: the weakbvarm and uniformirbvarm reduced-form
  priors, the Minnesota family via minnesotabvarm with glp marginal-likelihood
  tuning and hyperprior, and the structural apply functions svar.irf, svar.fevd,
  and historicalDecomposition. Use this skill WHENEVER the user mentions a
  Bayesian VAR, BVAR, SVAR, Structural VAR, a Minnesota/Litterman prior, a
  Normal-Inverse-Wishart or (semi)conjugate VAR prior, BVAR forecasting, or
  BVAR/SVAR impulse responses / FEVD / historical decomposition in MATLAB — even
  if they don't name the exact function. Also trigger when a user has
  macro/finance time series and asks for a shrinkage VAR, multivariate Bayesian
  forecasting, or structural shock analysis in MATLAB. Produces a runnable .m
  script covering the pipeline: prior → estimate → simulate → identify → apply.
metadata:
  author: Eduard Benet Cerda
  version: "3.0"
---

# Bayesian / Structural VAR (BVAR / SVAR) in MATLAB

This toolbox extends the Econometrics Toolbox so that **every** BVAR/SVAR
analysis is the same five stages:

1. **Choose a prior** — pick a prior model object.
2. **Estimate** — condition it on the data, getting a posterior object.
3. **Simulate** — draw reduced-form parameters from the posterior.
4. **Identify** — turn a covariance draw into a structural `impact` matrix.
5. **Apply** — IRF, FEVD, historical decomposition, forecasts.

The deliverable is a clean, runnable, commented `.m` file sectioned with `%%`
cells — not a chat explanation.

The unifying fact: **every prior object, built-in or from this toolbox, exposes
the same `estimate` / `simulate` / `forecast` surface.** Once stage 1 picks an
object, stages 2–3 are written the *same way* regardless of which prior it was.
Stages 4–5 operate on a `varm` object plus a structural `impact` matrix.

## Core rule: orchestrate, don't reimplement

The overriding rule for every line of code the skill produces.

- **Reuse in strict priority order: (1) an Econometrics Toolbox function → (2) a
  function in this toolbox (`tbx/svar`) → (3) only then write something new.** If
  a function exists, call it. Do **not** hand-roll a draw-to-`varm` converter,
  re-derive `Mu`/`V`, or write a manual IRF loop — use `minnesotabvarm`,
  `svar.varmFromCoefficients`, `svar.irf`. Almost every step below is one call.
- **Write sleek code; minimise `reshape` / `permute` / `squeeze`.** `svar.irf`
  and `svar.fevd` already return `(horizon+1) × series × shock`, so there is no
  `permute`/`reshape` dance — quantile or plot their output directly.
- **If something is genuinely missing**, a local function is acceptable, but give
  the user a one-paragraph **"should this live in the toolbox?"** note (reusable
  & general → propose adding it; one-off glue → keep it local). Extending the
  toolbox itself — new prior classes, new identification schemes, new apply
  functions — is a **development** task governed by the repo's `AGENTS.md`, not
  this skill.

---

## Stage 1 — Choose the prior

Pick the model object here. Everything downstream is identical across choices.

### Tier A — built-in Econometrics Toolbox priors (first choice)

Prefer the **direct constructors**; they take the full hyperparameter
specification as writable properties.

| Prior | Constructor | `estimate` returns | Set |
|---|---|---|---|
| Diffuse / noninformative | `diffusebvarm(m,p)` | `conjugatebvarm` (analytic) | nothing |
| Matrix-Normal-IW conjugate | `conjugatebvarm(m,p,…)` | `conjugatebvarm` (analytic) | `Mu,V,Omega,DoF` |
| Independent Normal-IW (semiconjugate) | `semiconjugatebvarm(m,p,…)` | `empiricalbvarm` (Gibbs) | `Mu,V,Omega,DoF` |
| Normal coeffs, fixed Σ | `normalbvarm(m,p,…)` | `normalbvarm` (analytic) | `Mu,V,Sigma` |

`m = numseries`, `p = numlags`. All accept `SeriesNames`, `IncludeConstant`,
`IncludeTrend`, `NumPredictors`. Do **not** use `bayesvarm` unless the user asks
for its Minnesota shortcut — it only sets these same properties.

### Tier B — this toolbox's priors

Ready-made, no hand-building. **Reduced-form improper priors** (thin
`conjugatebvarm` subclasses, data-free — pass `Y` at `estimate`):

```matlab
prior = weakbvarm(m, p);                              % Uhlig (2005) weak improper NIW
prior = uniformirbvarm(m, p);                         % Arias-Rubio-Ramirez-Waggoner uniform-IR
prior = uniformirbvarm(m, p, DeterminantShift=-1);
```

**Minnesota family** — the `minnesotabvarm` front door takes the data `Y` up
front and returns a complete prior (residual-variance scale and any dummy
observations already resolved, and the sample stored on the object):

```matlab
prior = minnesotabvarm(m, p, Y);                        % conjugate MNIW (default)
prior = minnesotabvarm(m, p, Y, lambda4=1, lambda5=1);  % + sum-of-coeff & dummy-initial-obs
prior = minnesotabvarm(m, p, Y, Type="inw", lambda2=0.5);
```

`Type` selects the family (aliases ignore hyphens/underscores/spaces):

| `Type` | Class | Posterior | Family-specific |
|---|---|---|---|
| `"mniw"` (default) | conjugate MNIW | analytic NIW | `lambda4`, `lambda5` |
| `"inw"` | independent Normal-IW | Gibbs | `lambda2` |
| `"normal"` | normal coeffs, fixed Σ | analytic | `lambda2` |

Common hyperparameters: `lambda1` overall tightness, `lambda3` lag decay,
`Vc` deterministic-coefficient variance, `PriorMean` (use `0` for differenced
series). `Psi` sets the residual-variance scale — a 1-by-`m` vector, or
`"exact"` (default) / `"conditional"` naming the estimator applied to `Y`. The
family's own constructor validates options, so passing `lambda2` to `"mniw"`
(where conjugacy pins it to 1) is an error, not a silent no-op. Full
hyperparameter detail: `references/priors.md`.

**Tune the tightnesses by marginal likelihood** (conjugate MNIW only) with
`glp` (Giannone–Lenza–Primiceri). A hyperparameter is *free* if you pass a
`hyperprior` or `[lo hi]` bounds, *fixed* if you pass a scalar:

```matlab
psi = svar.estimateResidualVariances(Y, p, Method="conditional");   % 1-by-m scale
[prior, info] = glp(m, p, Y, ...
    lambda1 = hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5]), ...   % free, Gamma prior
    lambda4 = [1e-4 50], ...                                        % free, flat prior
    lambda5 = 1, ...                                                % fixed
    Psi     = psi);
```

`glp` returns the prior built at the marginal-likelihood maximiser. For the full
hierarchical treatment — integrating *over* the hyperparameters rather than
fixing them at the mode — request a third output and a Metropolis chain:

```matlab
[prior, info, chain] = glp(m, p, Y, lambda1=..., Psi=psi, ...
    NumDraws=10000, BurnIn=10000, ProposalScale=0.8);
% chain.Coefficients / chain.Sigma already integrate over the hyperparameters
```

Feed `chain.Coefficients(:,:,d)` and `chain.Sigma(:,:,d)` straight into stage 5;
they skip stages 2–3. See `references/pipeline.md` for the hierarchical loop.

### Extending Tier B

If no built-in and no toolbox prior fits, building a **new** prior class is a
toolbox-development task — follow the recipes in the repo's `AGENTS.md`
(subclassing, hand-built `Mu`/`V`, Sims–Zha dummies). Do not inline a bespoke
prior class into an analysis script.

---

## Stage 2 — Estimate

Condition the prior on the data. Identical call for every prior object:

```matlab
[Posterior, Summary] = estimate(prior, Y, Display="off");
```

**Object-state rule.** Priors that already carry the sample — `minnesotabvarm`
and the `glp` output — are estimated with **no** `Y` argument (`estimate(prior)`);
passing a *different* `Y` errors, because their scale and dummies derive from the
stored sample. Data-free priors — every built-in, plus `weakbvarm` /
`uniformirbvarm` — **require** `Y`. The posterior class depends on the family:

| Prior | Posterior class | Draws come from |
|---|---|---|
| `conjugatebvarm`, `diffusebvarm`, `weakbvarm`, `uniformirbvarm`, Minnesota `mniw` | `conjugatebvarm` | analytic NIW |
| `normalbvarm`, Minnesota `normal` | `normalbvarm` | analytic (Σ fixed) |
| `semiconjugatebvarm`, Minnesota `inw` | `empiricalbvarm` | Gibbs sampler |

Set `rng(...)` before estimating any **sampler-based** (`inw` / semiconjugate)
model. Optional `X=` (exogenous predictors) and `Y0=` (presample) are accepted.

---

## Stage 3 — Simulate

Draw reduced-form parameters for structural work. Same call shape everywhere:

```matlab
rng(1);
[CoeffDraws, SigmaDraws] = simulate(Posterior, NumDraws=1000);
% CoeffDraws: numcoeff-by-NumDraws ;  SigmaDraws: m-by-m-by-NumDraws
```

The same **object-state rule** applies: call `simulate(Posterior, NumDraws=n)`
when the object embeds the sample (a posterior from `estimate`, a `minnesotabvarm`,
the `glp` output); call `simulate(prior, Y, NumDraws=n)` only when the object is
a data-free prior that still needs conditioning. Never pass `Y` to something that
already carries it — see `references/priors.md` §object-state.

For semiconjugate / `inw` models `estimate` returns an `empiricalbvarm` summary;
when you need raw draws, prefer `simulate(prior, Y, NumDraws=n)`.

---

## Stage 4 — Identify

Identification turns each reduced-form covariance draw `Sigma` into a structural
**impact matrix** `impact` (`numSeries × numShocks`) — the contemporaneous map
from structural shocks to responses.

The toolbox currently ships **recursive (Cholesky)** identification, which is one
line:

```matlab
impact = chol(Sigma, "lower");    % recursive ID; the SeriesNames order IS the assumption
```

**Always comment the variable ordering** wherever Cholesky is used — it *is* the
identifying assumption. To flip a shock's sign, negate that shock's `impact`
column (VARs are linear).

> Sign- and zero-restriction identification (Haar rotations, Arias–Rubio-Ramírez–
> Waggoner) is on the roadmap and will land as a dedicated identification API.
> Until then, supply `impact` yourself as above; a sign/zero routine you write
> for one analysis is local glue, and adding a scheme to the toolbox is an
> `AGENTS.md` task.

---

## Stage 5 — Apply

All apply functions take a `varm` object plus the `impact` matrix. Bridge from a
BVAR draw or posterior to a `varm` with the toolbox converters — never assemble a
`varm` by hand:

```matlab
Mdl = svar.varmFromCoefficients(Posterior, reshape(CoeffDraws(:,d), [], m), SigmaDraws(:,:,d)); % per draw
Mdl = bvar2var(Posterior);                                                      % whole model
```

`simulate` returns each coefficient draw as a `numcoeff×1` vector, so reshape it
to `m·k × m` (`reshape(CoeffDraws(:,d), [], m)`) before `svar.varmFromCoefficients`.
A `glp` chain already stores `chain.Coefficients(:,:,d)` in matrix form — pass it
straight through.

`svar.companionPower` is the expensive, impact-independent piece — compute it
**once per `varm`** and reuse it across `svar.irf` / `svar.fevd`:

```matlab
H   = 20;
Phi = svar.companionPower(Mdl, H);              % reduced-form MA blocks, reused below
ir  = svar.irf(Mdl,  impact, H, Phi=Phi);       % (H+1) × series × shock
fe  = svar.fevd(Mdl, impact, H, Phi=Phi);       % (H+1) × series × shock shares
```

`svar.irf` / `svar.fevd` already return `horizon × response × shock`, so quantile
or plot them directly — **no `permute`/`reshape`**.

`historicalDecomposition` takes a **BVAR object or a `varm`**, a **square**
`impact`, and `Y`, and returns a struct (`Contributions` `T−p × series × shock`,
`Total`, `StructuralShocks`, …). Passing the BVAR posterior directly is simplest:

```matlab
impact = chol(bvar2var(Posterior).Covariance, "lower");   % comment the ordering
HD     = historicalDecomposition(Posterior, impact, Y);   % or a varm from bvar2var
```

The canonical
credible-band loop (simulate → per-draw `svar.varmFromCoefficients` → identify →
`svar.irf`/`svar.fevd` → quantiles) and the FEVD/HD plotting patterns are in
`references/pipeline.md`.

**Forecasts** use the built-in surface directly on the posterior object:
`forecast(Posterior, fh, Y)` for mean/variance bands, or `simsmooth` with
NaN-padded data for full predictive paths (see `references/pipeline.md`).
Conditional forecasts are on the roadmap; until then, condition via `simsmooth`
with the known future values pinned in the NaN-padded sample.

---

## SSVS — what the toolbox actually offers

SSVS is **not** available for the BVAR (`*bvarm`) objects, built-in or toolbox.
It *is* available for Bayesian linear regression via `mixsemiconjugateblm` (and
custom priors via `customblm`). For a VAR: either approximate with heavy
shrinkage in a semiconjugate / `inw` Minnesota prior, or estimate
equation-by-equation with `mixsemiconjugateblm` and assemble the system. State
which you're doing; don't imply the BVAR objects support SSVS.

## Reference files
- `references/priors.md` — Stage 1 detail: the toolbox prior constructors, the
  `minnesotabvarm` / `glp` / `hyperprior` tuning API, the hyperparameter table,
  and the object-state rule. Read before choosing or tuning a prior.
- `references/pipeline.md` — Stages 2–5 detail: lag selection, forecast intervals
  / fan charts, the per-draw and hierarchical IRF/FEVD credible-band loops, FEVD
  plotting, and historical decomposition.

## Output conventions
- Orchestrate existing functions (Econ Toolbox → `tbx/svar` → new); no
  reimplementation of what ships. Any new local function comes with a
  "should this live in `tbx/svar`?" note; building new toolbox classes is an
  `AGENTS.md` task.
- Built-in prior first, then a toolbox prior.
- Set `rng` before any sampler-based (`inw` / semiconjugate) estimation or
  simulation, and before any `glp` Metropolis chain.
- Respect the object-state rule: pass `Y` only to data-free priors.
- Comment the variable ordering wherever Cholesky identification is used.
- Prefer natural output shapes; avoid `reshape`/`permute`/`squeeze` unless a
  function's documented output genuinely requires it.
- Section scripts with `%%` cell dividers so they run cell-by-cell.
