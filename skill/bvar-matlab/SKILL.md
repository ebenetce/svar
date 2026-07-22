---
name: bvar-matlab
description: >-
  Build, estimate, simulate, identify, and apply Bayesian and Structural Vector
  Autoregression (BVAR / SVAR) models in MATLAB. Covers the Econometrics Toolbox
  VAR prior objects (conjugatebvarm, semiconjugatebvarm, normalbvarm,
  diffusebvarm) AND this project's svar toolbox (tbx/svar): the weakbvarm and
  uniformirbvarm reduced-form priors, the Minnesota prior family via
  minnesotaSpec (with glp marginal-likelihood tuning and hyperprior), and the
  structural helpers svar.irf, svar.fevd, and historicalDecomposition. Use this
  skill WHENEVER the user mentions a Bayesian VAR, BVAR, SVAR, Structural VAR, a
  Minnesota/Litterman prior, a Normal-Inverse-Wishart or (semi)conjugate VAR
  prior, BVAR forecasting, or BVAR/SVAR impulse responses / FEVD / historical
  decomposition in MATLAB — even if they don't name the exact function. Also
  trigger when a user has macro/finance time series and asks for a shrinkage VAR,
  multivariate Bayesian forecasting, or structural shock analysis in MATLAB.
  Produces a runnable .m script covering the full pipeline: prior → estimate →
  simulate → identify → apply.
metadata:
  author: Eduard Benet Cerda
  version: "2.0"
---

# Bayesian / Structural VAR (BVAR / SVAR) in MATLAB

This skill builds the **full BVAR/SVAR pipeline** as five explicit stages:

1. **Choose prior** — pick a model object.
2. **Estimate** — condition on the data.
3. **Simulate** — draw reduced-form parameters.
4. **Identify** — turn a covariance draw into a structural `impact` matrix.
5. **Apply** — IRF, FEVD, historical decomposition.

The deliverable is a clean, runnable, commented `.m` file sectioned with `%%`
cells — not a chat explanation.

The unifying fact: **every prior object — built-in or this project's — inherits
the same `estimate` / `simulate` / `forecast` surface.** So once stage 1 picks an
object, stages 2–3 are written the *same way* regardless of which prior it is.
Stages 4–5 operate on a `varm` object plus a structural `impact` matrix.

## Core rule: orchestrate, don't reimplement

This is the overriding rule for every line of code the skill produces.

- **Reuse in strict priority order: (1) an Econometrics Toolbox function → (2) a
  function in this project (`tbx/svar`) → (3) only then write something new.** If
  a function already exists on that path, call it. Do **not** re-derive `Mu`/`V`,
  hand-roll a `bvarDrawToVARM`, or write a manual IRF loop — use
  `minnesotaSpec`, `varmFromCoefficients`, `svar.irf`. Almost every step below is
  one function call.
- **Write sleek code; minimize `reshape` / `permute` / `squeeze`.** Prefer the
  natural output shapes the existing functions already return. Concretely:
  `svar.irf` and `svar.fevd` already return `(horizon+1) × series × shock`, so
  there is no `permute([1 3 2])` or `reshape` dance to write — the built-in
  `varm` workflow needs those; this toolbox does not.
- **If something is genuinely missing** (nothing in Econ Toolbox and nothing in
  `tbx/svar` fits), a local function is acceptable — but (a) write it in a shape
  fit for `tbx/svar` (same arg conventions, `varm`/`bvar` interfaces, package
  placement), and (b) give the user a one-paragraph **"should this live in the
  toolbox?" assessment**: yes/no with reasons (reusable & general → add it to
  `tbx/svar`; one-off glue for this analysis → keep it local). Never silently
  accrete throwaway helpers.

---

## Stage 1 — Choose the prior

Pick the model object here. Everything downstream is identical across choices.

### Tier A — built-in Econometrics Toolbox priors (first choice)

Prefer the **direct constructors**, which take the full hyperparameter
specification as writable properties:

| Prior | Constructor | `estimate` returns | Set |
|---|---|---|---|
| Diffuse / noninformative | `diffusebvarm(m,p)` | `conjugatebvarm` (analytic) | nothing |
| Matrix-Normal-IW conjugate | `conjugatebvarm(m,p,…)` | `conjugatebvarm` (analytic) | `Mu,V,Omega,DoF` |
| Independent Normal-IW (semiconjugate) | `semiconjugatebvarm(m,p,…)` | `empiricalbvarm` (Gibbs) | `Mu,V,Omega,DoF` |
| Normal coeffs, fixed Σ | `normalbvarm(m,p,…)` | `normalbvarm` (analytic) | `Mu,V,Sigma` |

`m = numseries`, `p = numlags`. All accept `SeriesNames`, `IncludeConstant`,
`IncludeTrend`, `NumPredictors`.

**Do not use `bayesvarm`** unless the user explicitly asks for its Minnesota
shortcut syntax — it is only a convenience wrapper that sets the same properties.

### Tier B — this project's priors (`tbx/svar/priors/`)

When the user wants a Minnesota/Litterman prior or a research-grade reduced-form
prior, these are ready-made — no hand-building required.

**Reduced-form improper priors** (thin `conjugatebvarm` subclasses):
```matlab
prior = weakbvarm(m, p);                 % Uhlig (2005) weak improper NIW
prior = uniformirbvarm(m, p);            % Arias–Rubio-Ramírez–Waggoner uniform-IR
prior = uniformirbvarm(m, p, DeterminantShift=-1);
```
Both accept the standard `SeriesNames`/`IncludeConstant`/`IncludeTrend`/
`NumPredictors` options.

**Minnesota family** — build via the `minnesotaSpec` factory, which needs a
residual-variance scale `psi` from `estimateResidualVariances`:
```matlab
psi   = estimateResidualVariances(Y, p);          % 1-by-m variance scale
prior = minnesotaSpec("mniw").build(m, p, psi);   % conjugate MNIW (default)
% Methods: "mniw"|"conjugate", "inw"|"semiconjugate", "normal"|"fixedsigma"
```
`minnesotaSpec` is the public entry point; the concrete `minnesota*bvarm` classes
are hidden. Set hyperparameters on the spec (`lambda1` tightness, `lambda3` lag
decay, `lambda4`/`lambda5` sum-of-coefficients & dummy-initial-obs for `mniw`,
`lambda2` cross-tightness for `inw`/`normal`, `Vc` deterministic variance). See
`references/priors.md` for the full hyperparameter table.

**Tune the hyperparameters by marginal likelihood** (conjugate MNIW only) with
`glp` (Giannone–Lenza–Primiceri); make a hyperparameter "free" by passing a
`hyperprior` or `[lo hi]` bounds instead of a scalar:
```matlab
[prior, info] = glp(m, p, Y, psi);        % tune lambdas by marginal likelihood
spec = minnesotaSpec("mniw", lambda1=hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5]));
[prior, info] = glp(m, p, Y, psi, Spec=spec);
```

### Tier C — nothing fits (last resort)

Only if no built-in and no `tbx/svar` prior matches: subclass the closest
built-in and override the constructor to reset `Mu`/`V`/`Omega`/`DoF` (or
`Sigma`), or hand-build `Mu`/`V`. This is exactly how the Tier-B priors were
built — read them first. The recipes (λ layout, Minnesota `Mu`/`V`, Sims–Zha
dummies, what to override per prior type) are in `references/priors.md`. If you
write a new prior class, give the "should this live in `tbx/svar`?" assessment.

---

## Stage 2 — Estimate

Condition the prior on the data. Identical call for every prior object:
```matlab
[Posterior, Summary] = estimate(prior, Y, Display="off");
```
`estimate` returns a **posterior model object**, whose class depends on the
family:

| Prior | Posterior class | Draws come from |
|---|---|---|
| `conjugatebvarm`, `diffusebvarm`, `weakbvarm`, `uniformirbvarm`, Minnesota `mniw` | `conjugatebvarm` | analytic NIW |
| `normalbvarm`, Minnesota `normal` | `normalbvarm` | analytic (Σ fixed) |
| `semiconjugatebvarm`, Minnesota `inw` | `empiricalbvarm` | Gibbs sampler |

Set `rng(...)` before estimating any **sampler-based** (semiconjugate / `inw`)
model. Optional `X=` (exogenous predictors), `Y0=` (presample) are accepted.

---

## Stage 3 — Simulate

Draw reduced-form parameters for structural work. Same call shape everywhere:
```matlab
rng(1);
[CoeffDraws, SigmaDraws] = simulate(prior, Y, NumDraws=1000);
% CoeffDraws: numcoeff-by-NumDraws ;  SigmaDraws: m-by-m-by-NumDraws
```

**Object-state rule for `simulate`:** pass `Y` when the object is still a
**prior** that must be conditioned on the sample. Omit `Y` and call
`simulate(Posterior, NumDraws=n)` only when the object already **is** a posterior
that embeds the sample. Never pass `Y` to something that already carries it — see
`references/priors.md` §8.

For semiconjugate / `inw` models, `estimate` returns an `empiricalbvarm` summary;
when you need raw draws for stage 5, prefer `simulate(prior, Y, NumDraws=n)`.

---

## Stage 4 — Identify

Identification turns each reduced-form covariance draw `Sigma` into a structural
**impact matrix** `impact` (`numSeries × numShocks`), the mapping from structural
shocks to contemporaneous responses.

There is **no `identify` function yet** — supply `impact` yourself. The recursive
(Cholesky) scheme in `SeriesNames` order is one line:
```matlab
impact = chol(Sigma, "lower");    % recursive identification; comment the ordering
```
Always comment the variable ordering wherever Cholesky is used — it defines the
identifying assumption.

> A dedicated `svar.Identification.identify` (sign restrictions, Haar rotations,
> etc.) is planned; when it lands it supersedes the manual `impact` construction.

---

## Stage 5 — Apply

All apply functions take a `varm` object plus the `impact` matrix. Bridge from a
BVAR draw or posterior to a `varm` with the project's converters — do not
assemble a `varm` by hand:

```matlab
Mdl = varmFromCoefficients(prior, CoeffDraws(:,d), SigmaDraws(:,:,d)); % per draw
Mdl = bvar2var(Posterior);                                            % whole model
```

Then compute structural objects. `svar.companionPower` is the expensive,
impact-independent piece — compute it **once per `varm`** and reuse it:
```matlab
H   = 20;
Phi = svar.companionPower(Mdl, H);              % reduced-form MA blocks, reused below
ir  = svar.irf(Mdl,  impact, H, Phi=Phi);       % (H+1) × series × shock, directly
fe  = svar.fevd(Mdl, impact, H, Phi=Phi);       % (H+1) × series × shock shares
HD  = historicalDecomposition(Mdl, impact, Y);  % impact must be square (m×m)
```

`svar.irf` / `svar.fevd` already return `horizon × response × shock`, so plot or
quantile them directly — **no `permute`/`reshape`**. The canonical draw loop for
posterior credible bands (simulate → per-draw `varmFromCoefficients` → identify →
`svar.irf`/`svar.fevd` → quantiles), forecast fan charts, and the
`historicalDecomposition` struct fields are in `references/pipeline.md`.

To flip a shock's sign, negate that shock's `impact` column (VARs are linear).

---

## SSVS — what the toolbox actually offers

SSVS is **not** available for the BVAR (`*bvarm`) objects, built-in or project.
It *is* available for Bayesian linear regression via `mixsemiconjugateblm` (and
arbitrary custom priors via `customblm`). For a VAR: either approximate with
heavy shrinkage in a semiconjugate / `inw` Minnesota prior, or estimate
equation-by-equation with `mixsemiconjugateblm` and assemble the system. State
which you're doing; don't imply the BVAR objects support SSVS.

## Reference files
- `references/priors.md` — Stage 1 detail: project-prior constructors and the
  `minnesotaSpec`/`glp`/`hyperprior` tuning API, the hyperparameter table, and
  (demoted, last-resort) the λ/`Mu`/`V` layout, hand-built Minnesota moments,
  Sims–Zha dummies, subclassing recipes, and the object-state / diffuse-OLS-NIW
  rules. Read before choosing or building a prior.
- `references/pipeline.md` — Stages 2–5 detail: lag selection, forecast
  intervals / fan charts, the per-draw IRF/FEVD credible-band loop, FEVD
  plotting, and historical decomposition.

## Output conventions
- Orchestrate existing functions (Econ Toolbox → `tbx/svar` → new); no
  reimplementation of what ships. Any new local function comes with a
  "should this live in `tbx/svar`?" assessment.
- Built-in prior first, then a `tbx/svar` prior; subclass/hand-build only as a
  last resort.
- Set `rng` before any sampler-based (`inw` / semiconjugate) estimation or
  simulation.
- Comment the variable ordering wherever Cholesky identification is used.
- Prefer natural output shapes; avoid `reshape`/`permute`/`squeeze` unless a
  function's documented output genuinely requires it.
- Section scripts with `%% ` cell dividers so they run cell-by-cell.
