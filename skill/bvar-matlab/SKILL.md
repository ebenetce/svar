---
name: bvar-matlab
description: >-
  Build, estimate, forecast, and structurally analyze Bayesian Vector
  Autoregression (BVAR) models in MATLAB using the Econometrics Toolbox VAR
  prior objects (conjugatebvarm, semiconjugatebvarm, normalbvarm, diffusebvarm)
  and their subclasses. Use this skill WHENEVER the user mentions a Bayesian
  VAR, BVAR, SVAR, Structural VAR, a Minnesota/Litterman prior, a Normal-Inverse-Wishart or
  (semi)conjugate VAR prior, a custom or extended BVAR prior, BVAR forecasting,
  or BVAR impulse responses / FEVD / historical decomposition in MATLAB — even
  if they don't name the exact function. Also trigger when a user has
  macro/finance time series and asks for a shrinkage VAR, multivariate Bayesian
  forecasting, or structural shock analysis in MATLAB. Produces a runnable .m
  script (or class file) covering the full pipeline: prior to estimation to
  forecast to IRF/FEVD.
metadata:
  author: Eduard Benet Cerda
  version: "1.0"
---

# Bayesian VAR (BVAR) in MATLAB — Econometrics Toolbox

This skill builds the **full BVAR pipeline** in MATLAB's Econometrics Toolbox:
prior specification → posterior estimation → forecasting → structural analysis
(IRF, FEVD). The deliverable is a clean, runnable, commented `.m` file (script,
function, or class), not a chat explanation.

## Core design rule: use the direct constructors, not `bayesvarm`

**Prefer the direct prior-object constructors** — `conjugatebvarm`,
`semiconjugatebvarm`, `normalbvarm`, `diffusebvarm`. Each one accepts the
**full hyperparameter specification** (`Mu`, `V`, `Omega`, `DoF`, `Sigma`) as
writable properties, so you get complete control over the prior.

`bayesvarm` is **only** a convenience wrapper that sets those same properties
from a handful of Minnesota shortcut args (`Center/SelfLag/CrossLag/Decay`). It
is not flexible enough for real work: it cannot express an arbitrary `Mu`/`V`,
asymmetric own-vs-cross treatment, or block structures. **Do not use
`bayesvarm`** unless the user explicitly asks for the Minnesota-shortcut
convenience syntax. Build `Mu` and `V` yourself and pass them to the direct
constructor instead.

## Choosing the constructor

| Prior | Constructor | Posterior from `estimate` | Key hyperparameters you set |
|---|---|---|---|
| Diffuse / noninformative | `diffusebvarm(m,p)` | analytic `conjugatebvarm` | none (flat) |
| Matrix-Normal-Inverse-Wishart conjugate | `conjugatebvarm(m,p,...)` | analytic `conjugatebvarm` | `Mu`, `V`, `Omega`, `DoF` |
| Independent Normal-Inverse-Wishart (semiconjugate) | `semiconjugatebvarm(m,p,...)` | **Gibbs** → `empiricalbvarm` | `Mu`, `V`, `Omega`, `DoF` |
| Normal coefficients, fixed Σ | `normalbvarm(m,p,...)` | `normalbvarm` (analytic given Σ) | `Mu`, `V`, `Sigma` |

`m = numseries`, `p = numlags`. All four also accept `SeriesNames`,
`IncludeConstant`, `IncludeTrend`, `NumPredictors`.

The coefficient vector is `λ = vec([Φ1 Φ2 … Φp c δ Β]')`. `Mu` is the prior
mean of `λ`; `V` is its prior covariance factor. Build them explicitly — see
`references/priors.md` for the layout of `λ` and how to place Minnesota-style
shrinkage into `Mu`/`V` by hand (which is what `bayesvarm` does internally, but
under your control).

## Extending priors: subclass first, new class only as a last resort

**Decision rule when no built-in prior is an exact match for what the user
wants:**

1. **Exact built-in match** → use the direct constructor (above). Done.
2. **No exact match (the default case for any non-standard prior)** →
   **immediately subclass the closest built-in object and override only what
   differs.** This is the *first* thing to reach for, not a fallback. You
   inherit the sampler, the data handling, the forecast recursion, and all the
   IRF/FEVD/decomposition plumbing for free, and change only the piece that is
   actually different.
3. **Fully custom class from scratch (no `*bvarm` superclass)** → reach for
   this **only if** subclassing genuinely cannot express the prior — e.g. the
   likelihood itself changes, or the parameter block does not map onto the NIW
   structure at all. When you make this call, say so explicitly and explain why
   subclassing was insufficient.

In practice almost every "the toolbox doesn't ship this prior" request is
solved at step 2. Do **not** jump to a brand-new class because the prior looks
unusual — try subclassing first.

### What to override when subclassing

These objects are MATLAB classes. Pick the closest built-in as the superclass
and override the minimum:

- **Just the prior moments** (asymmetric conjugate prior, custom shrinkage,
  block structure, hierarchical hyperprior defaults) → override the
  **constructor** to reset `Mu` / `V` / `Omega` / `DoF` (or `Sigma`) after
  delegating to the superclass. Nothing else needs to change.
- **Dummy-observation priors** (sum-of-coefficients, Sims–Zha, single-unit-root)
  → override `estimate` / `simulate` to inject the dummy observations before
  delegating to the superclass.
- **A different sampler or posterior mechanics** → override `estimate` (and
  `simulate` if draws are needed), then let `forecast` / `irf` / `fevd` operate
  on the resulting object.
- **Changed forecast or structural recursion** → override `forecast` /
  `simulate` only.

Pattern:
```matlab
classdef myPriorBVARM < semiconjugatebvarm
    % Custom BVAR prior extending the semiconjugate object.
    methods
        function obj = myPriorBVARM(numseries, numlags, varargin)
            obj = obj@semiconjugatebvarm(numseries, numlags, varargin{:});
            % override Mu / V / Omega / DoF here to encode the custom prior
        end
        % Override estimate/forecast/simulate ONLY for the piece that differs
        % (e.g. inject dummy observations, then call the superclass method).
    end
end
```
Guidance on what to override for each extension (and the `customblm`-style
fully-custom approach for the rare step-3 case) is in `references/priors.md`.
Deliver a `classdef` file subclassing the closest built-in.

## Workflow

### 1. Clarify only what you must
If the user gave you dimensions, lags, and a prior, proceed. If the prior is
ambiguous ("a good forecasting BVAR"), default to `semiconjugatebvarm` with
Minnesota-style `Mu`/`V` (own-lag-1 mean = 1 for levels / 0 for differenced
data) and say so in a comment. Don't over-interview.

### 2. Build the prior with a direct constructor
```matlab
m = size(Y,2); p = 4;
names = ["INFL" "UNRATE" "FEDFUNDS"];

PriorMdl = semiconjugatebvarm(m, p, SeriesNames=names);

% Set the full prior explicitly (see references/priors.md for Mu/V layout):
PriorMdl.Mu    = buildMinnesotaMu(m, p);          % prior coeff mean
PriorMdl.V     = buildMinnesotaV(m, p, lambda);   % prior coeff covariance
PriorMdl.Omega = diag(resVar);                    % IW scale
PriorMdl.DoF   = m + 2;                           % IW dof (finite mean)
```
(`buildMinnesotaMu`/`buildMinnesotaV` are helper functions you write; the
construction is spelled out in `references/priors.md`.)

If the prior is non-standard, build it as a subclass per the decision rule
above rather than contorting the built-in shortcuts.

### 3. Estimate or condition on the data
```matlab
rng(1);  % required for semiconjugate/normal (sampler-based)
PosteriorMdl = estimate(PriorMdl, Y, Display="off");
summarize(PosteriorMdl);
```
- `diffuse`/`conjugate`: analytic posterior, returns a `conjugatebvarm`.
- `semiconjugate`: Gibbs sampler → `empiricalbvarm` (draws). Tune with
  `NumDraws`, `BurnIn`, `Thin`.
- `normal`: analytic given fixed Σ, returns a `conjugatebvarm`.

For data-conditioned parameter draws, call `simulate` on the prior object and
pass the sample:
```matlab
[CoeffDraws, SigmaDraws] = simulate(PriorMdl, Y, NumDraws=1000);
```
This is the right path for semiconjugate BVARs when the next step needs
coefficient/covariance draws for structural analysis. `estimate` is useful for
summaries, but its semiconjugate output is an `empiricalbvarm`, which is not a
drop-in replacement for a `varm` in IRF/FEVD work.

**Reduced-form SVAR draws? Use the object-state rule.** When the goal is
reduced-form parameter draws for external identification (sign restrictions,
Haar rotations), first decide whether the object is still an unconditioned prior
or already represents a posterior. For any built-in prior (`diffusebvarm`,
`conjugatebvarm`, `semiconjugatebvarm`, `normalbvarm`) or subclass, call
`simulate(PriorMdl, Y, NumDraws=n)` when `PriorMdl` is a pure prior that still
needs to be conditioned on `Y`; call `simulate(PosteriorMdl, NumDraws=n)` only
when `PosteriorMdl` is a posterior object that supports `simulate`. If a
posterior object already embeds the estimation data (for example a custom NIW
Minnesota posterior or a manually constructed diffuse-prior NIW posterior), do
**not** pass `Y` again. See `references/priors.md` §8 for the full decision rule
and the diffuse-prior OLS NIW shortcut.

### 4. Forecast
BVAR `forecast` conditions on the supplied estimation sample. Forecasting from
the prior object is valid and often preferable for semiconjugate models:
```matlab
fh = 12;
[YF, YFStd] = forecast(PriorMdl, fh, Y);       % YF: fh-by-m
```
Intervals / fan charts and path simulation: `references/pipeline.md`.

### 5. Structural analysis (IRF / FEVD)
Do **not** call `irf` or `fevd` directly on `*bvarm` objects. Draw or summarize
the BVAR coefficients/covariance, construct a compatible `varm` object, then
call the `varm` IRF/FEVD methods. Cholesky-orthogonalized results use the
`SeriesNames` order, so comment the ordering.
```matlab
numDraws = 1000;
[CoeffDraws, SigmaDraws] = simulate(PriorMdl, Y, NumDraws=numDraws);

H = 20;
for d = 1:numDraws
    Mdl_d = bvarDrawToVARM(PriorMdl, CoeffDraws(:,d), SigmaDraws(:,:,d));
    Response_d = irf(Mdl_d, NumObs=H);        % horizon-by-response-by-shock
    Decomp_d   = fevd(Mdl_d, NumObs=H);
end
```
The `bvarDrawToVARM` helper, coefficient-map unpacking, IRF credible bands, FEVD
plotting, and the historical-decomposition caveat are in
`references/pipeline.md`.

## SSVS — what the toolbox actually offers
SSVS is **not** available for the BVAR (`*bvarm`) objects. It *is* available for
Bayesian linear regression via `mixsemiconjugateblm` (and arbitrary custom
priors via `customblm`). So for a VAR you have two honest options: (a) approximate
with heavy shrinkage in `semiconjugatebvarm`, or (b) estimate the VAR
equation-by-equation with `mixsemiconjugateblm` and assemble the system
yourself. State which you're doing; don't imply the BVAR objects support SSVS.

## Reference files
- `references/priors.md` — the `λ` coefficient layout, how to build `Mu`/`V`
  for Minnesota and other shrinkage by hand, every hyperparameter, the Sims–Zha
  dummy-observation construction, and subclassing recipes (what to override per
  prior type). Read before building `Mu`/`V`, tuning, or subclassing.
- `references/pipeline.md` — lag selection, forecast intervals/fan charts,
  per-draw IRF credible bands, FEVD plotting, historical-decomposition caveat.

## Output conventions
- Use direct constructors; only use `bayesvarm` if the user asks for it.
- No exact built-in prior match → subclass the closest built-in and override
  only the methods that differ (constructor for `Mu`/`V`/`Omega`/`DoF`;
  `estimate`/`simulate`/`forecast` as needed). Write a fully custom class only
  when subclassing genuinely can't express the prior, and say why.
- Set `rng` before any sampler-based estimation.
- Comment the variable ordering wherever Cholesky identification is used.
- For custom priors, deliver a `classdef` file subclassing the closest built-in.
- Section scripts with `%% ` cell dividers so they run cell-by-cell.
