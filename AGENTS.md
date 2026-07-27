# AGENTS.md — developing the `svar` toolbox

Guidance for changing this repo's code. For *using* the toolbox to build a model,
the `bvar-matlab` skill (`.agents/skills/bvar-matlab/`) is the reference; this
file is about *extending* it.

## Mandatory MATLAB project setup

Before doing any MATLAB work in this repository, open the MATLAB project:

```matlab
openProject(fullfile(pwd,"SVAR.prj"));
```

If you are not starting from the repository root, open the project with an
absolute path to `SVAR.prj`.

Do not run tests, examples, analysis, or edits that depend on MATLAB project
configuration until the project is open. If the project cannot be opened, stop
and report that blocker instead of continuing with an unconfigured MATLAB
session.

## What this repo is

`svar` is an **extension of the MATLAB Econometrics Toolbox** for Bayesian and
Structural VARs. It does not replace the toolbox — every prior object subclasses
a built-in `*bvarm` class and keeps the built-in `estimate` / `simulate` /
`forecast` surface. The whole toolbox is organised around one pipeline, and new
code should slot into one of its stages:

1. **choose prior** → prior model objects (`tbx/svar/priors/`)
2. **estimate** → inherited from the built-in base classes
3. **simulate** → inherited
4. **identify** → produce a structural `impact` matrix (Cholesky today)
5. **apply** → `svar.irf`, `svar.fevd`, `historicalDecomposition`, forecasts

## Layout & packaging

| Path | Contents | Callable as |
|---|---|---|
| `tbx/svar/+svar/` | array-returning structural / numeric helpers | `svar.irf`, `svar.fevd`, `svar.companionPower`, `svar.companionMatrix`, `svar.varmFromCoefficients`, `svar.estimateResidualVariances`, `svar.logDetPD` |
| `tbx/svar/` | user-facing orchestration functions | `bvar2var`, `hyperprior`, `logMarginalLikelihood`, `marginalLikelihood`, `historicalDecomposition` |
| `tbx/svar/priors/` | prior **front doors** (on the path) | `minnesotabvarm`, `weakbvarm`, `uniformirbvarm`, `glp` |
| `tbx/svar/priors/+svar/` | hidden concrete Minnesota classes + abstract base | `svar.minnesotamniwbvarm`, `svar.minnesotainwbvarm`, `svar.minnesotanbvarm`, `svar.minnesotabvarmBase` |
| `tbx/doc/` | one `*.md` per public function/class | — |
| `tests/` | `matlab.unittest` classes, `*Test.m` | — |

**Packaging rule:** a numeric helper that returns arrays and is composed inside
other functions goes in `+svar`. A front door a user calls directly (a prior
constructor, a tuner) is bare on the path. Concrete prior classes are hidden in
`priors/+svar`; users reach them through the `minnesotabvarm` dispatcher.

## Class hierarchy

```
conjugatebvarm (built-in)
├── weakbvarm                    reduced-form improper NIW (Uhlig 2005)
├── uniformirbvarm               uniform on the impact rotation
└── svar.minnesotamniwbvarm      + svar.minnesotabvarmBase   (Type="mniw")
semiconjugatebvarm (built-in)
└── svar.minnesotainwbvarm       + svar.minnesotabvarmBase   (Type="inw")
normalbvarm (built-in)
└── svar.minnesotanbvarm         + svar.minnesotabvarmBase   (Type="normal")
```

`svar.minnesotabvarmBase` (abstract) holds the shared Minnesota state: the stored
sample `Y`, `ResidualVariances`, `lambda1`/`lambda3`/`Vc`/`PriorMean`, and the
`configureMinnesota` / `resolvePsi` / `buildMinnesotaPriorMean` helpers.

## Invariants (easy to get wrong)

- **Materialise the full prior in the constructor.** `minnesotabvarm` folds the
  sum-of-coefficients / dummy-initial-observation pseudo-observations into
  `Mu`/`V`/`Omega`/`DoF` at construction (`applyDummies` → `updateNIW`). That is
  why `logMarginalLikelihood` and `estimate` need **no** Minnesota special case.
  Any new prior that adds pseudo-observations must do the same — never patch the
  dummies in downstream, or every consumer has to remember to apply them.
- **DoF = m+2.** The Minnesota classes override the inherited `conjugatebvarm`
  DoF default; `V` is built for `DoF = m+2` (the (DoF−m−1) scaling is then 1).
  Changing one without the other silently rescales the prior.
- **Conjugacy pins `lambda2 = 1`.** The `Σ⊗V` Kronecker structure of the
  conjugate family cannot represent separate own/cross tightness. A free
  `lambda2` requires the `inw` (semiconjugate) or `normal` family. Do not add
  `lambda2` to `mniw`.
- **Sample storage.** Minnesota priors store `Y` and default `estimate`/
  `simulate`/`forecast` to it, rejecting a *foreign* sample. Reduced-form priors
  (`weakbvarm`, `uniformirbvarm`) are data-free. Preserve this split.
- **`glp` is conjugate-`mniw` only.** It maximises the *analytic* marginal
  likelihood; `inw`/`normal`/semiconjugate have none. `logMarginalLikelihood`
  errors on those on purpose — keep it that way.
- **Output shape `(H+1) × series × shock`.** All apply functions return this.
  `svar.companionPower(Mdl, H)` is the impact-independent MA piece; new apply
  functions accept `Phi=` so a caller can compute it once and reuse it.

## Extending — read the canonical example, then match it

Do not invent a new shape. Each extension type has a reference implementation;
read it and follow its conventions (arguments block, error IDs `func:reason`,
name-value forwarding via `namedargs2cell`).

- **New reduced-form prior** (improper NIW variant): subclass `conjugatebvarm`,
  set `Mu`/`V`/`Omega`/`DoF` in the constructor. Template: `weakbvarm.m` (minimal)
  or `uniformirbvarm.m` (with an option).
- **New Minnesota-family variant**: subclass the matching built-in **and**
  `svar.minnesotabvarmBase`; call `configureMinnesota` for the shared state, then
  build your moments. Template: `priors/+svar/minnesotamniwbvarm.m`. Register a
  `Type` alias in the `minnesotabvarm` front door.
- **New identification scheme** (sign / zero restrictions are next): the contract
  is a map from a reduced-form `Sigma` draw to an `impact` matrix
  (`numSeries × numShocks`); set-identified schemes return a set / accept a
  rotation. The apply functions already take `impact` and stay agnostic, so a new
  scheme is additive — put it under a `+svar` identification package and leave
  `svar.irf`/`svar.fevd` untouched.
- **New apply function**: signature `(varMdl, impact, horizon, nvp.Phi)`, return
  `(horizon+1) × series × shock`, reuse `svar.companionPower`. Template:
  `+svar/irf.m`, `+svar/fevd.m`.

### Minnesota `Mu` / `V` layout (stable reference)

Per-equation coefficient block is `[Φ1 … Φp c δ B]'`; the estimable vector is
`λ = vec(Λ')` stacking equations, `k = m·p + IncludeConstant + IncludeTrend +
NumPredictors` coefficients per equation, so `λ` has `m·k` elements.

- `Mu` — prior mean, `PriorMean` on the own first lag, 0 elsewhere.
- `V` — coefficient covariance factor (`Σ⊗V` conjugate; standalone for semiconj.).
  Lag block: `Var(lag ℓ, regressor j) = lambda1² / (ℓ^(2·lambda3) · ψ_j)`;
  deterministic/exogenous columns get `Vc`.
- `Omega = diag(ψ)`, `DoF = m+2`.

Always verify indexing against `PriorMdl.Mu` on a small example — the ordering is
construction-dependent.

### Recognising legacy NIW samplers

Hand-rolled reduced-form SVAR code usually draws Σ ~ IW then β | Σ ~ N. That is
exactly a `conjugatebvarm` with `Omega = T·Σ̂`, `DoF = T`, `Mu = OLS coeffs`,
`V = inv(X'X)`, plus `simulate`. Offer the object equivalent instead of porting
the loop — but confirm the coefficient ordering matches `PriorMdl.Mu`.

## Testing, docs, build

- **Test:** `buildtool test` (default), or `runtests("tests", IncludeSubfolders=true)`,
  or the MATLAB MCP `run_matlab_test_file`. Tests are `matlab.unittest.TestCase`
  classes named `<thing>Test.m`; a `ProjectFixture` loads the project. Add tests
  for every new public function (shape, invariants, error IDs).
- **Check:** `buildtool check` runs `CodeIssuesTask` on `tbx` with **zero**
  warning/info tolerance — new code must be clean. Also available via MCP
  `check_matlab_code` per file.
- **Docs:** one `tbx/doc/<name>.md` per public function/class, plus in-file help
  (H1 line, `See also`). `buildtool doc` converts the markdown via DocMaker.
  A new public function needs a `.md` and a `helptoc.md` entry.
- **Reuse order for any new code:** Econometrics Toolbox → existing `tbx/svar` →
  only then something new. Assess whether a new helper is general enough to live
  in `tbx/svar` (add it) or is one-off glue (keep it in the caller).
