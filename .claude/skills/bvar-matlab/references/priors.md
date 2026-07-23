# BVAR Priors & Hyperparameters — Usage Reference

Stage-1 detail for *choosing and configuring* a prior with the existing
constructors. Building a **new** prior class (hand-built `Mu`/`V`, subclassing,
Sims–Zha dummy construction, the NIW Rosetta stone) is a toolbox-development
task — see the repo `AGENTS.md`.

## Table of contents
1. Constructor cheat-sheet
2. The Minnesota front door (`minnesotabvarm`)
3. Hyperparameters
4. Choosing hyperparameter values
5. Tuning by marginal likelihood (`glp` + `hyperprior`)
6. The object-state rule (who carries the sample)
7. SSVS — what is and isn't available

---

## 1. Constructor cheat-sheet

| Constructor | Prior | `estimate` → | You set |
|---|---|---|---|
| `diffusebvarm(m,p)` | flat / noninformative | analytic `conjugatebvarm` | nothing |
| `conjugatebvarm(m,p,…)` | Matrix-Normal-IW (Σ⊗V) | analytic `conjugatebvarm` | `Mu,V,Omega,DoF` |
| `semiconjugatebvarm(m,p,…)` | independent Normal-IW | Gibbs → `empiricalbvarm` | `Mu,V,Omega,DoF` |
| `normalbvarm(m,p,…)` | normal coeffs, Σ fixed | analytic `normalbvarm` | `Mu,V,Sigma` |
| `weakbvarm(m,p)` | Uhlig (2005) weak improper NIW | analytic `conjugatebvarm` | nothing |
| `uniformirbvarm(m,p)` | uniform on the impact rotation | analytic `conjugatebvarm` | `DeterminantShift` |
| `minnesotabvarm(m,p,Y,…)` | Minnesota / Litterman family | depends on `Type` | `lambda*`, `Vc`, `Psi`, … |

The built-in four and `weakbvarm` / `uniformirbvarm` are **data-free**: you pass
`Y` at `estimate`. `minnesotabvarm` takes `Y` **up front** and stores it (see §6).

`uniformirbvarm(m, p, DeterminantShift=-1)` shifts the determinant term when
matching a specific reduced-form convention; the default is fine otherwise.

Conjugacy constraint for the conjugate family: the coefficient covariance is
`Σ⊗V`, so own- and cross-variable tightness cannot differ (that is why `mniw`
has no free `lambda2`). If you need asymmetric own-vs-cross shrinkage, use the
`"inw"` (semiconjugate) or `"normal"` family, which carry a free `lambda2`.

---

## 2. The Minnesota front door (`minnesotabvarm`)

```matlab
prior = minnesotabvarm(numseries, numlags, Y, Name=Value)
```

`Type` selects the family; aliases ignore hyphens, underscores, and spaces:

| `Type` | Aliases | Class | Posterior | Family-specific |
|---|---|---|---|---|
| `"mniw"` (default) | `conjugate`, `matrixnormal` | conjugate MNIW | analytic NIW | `lambda4`, `lambda5` |
| `"inw"` | `independent`, `kadiyala`, `semiconjugate` | independent Normal-IW | Gibbs | `lambda2` |
| `"normal"` | `litterman`, `fixed`, `fixedsigma` | normal coeffs, fixed Σ | analytic | `lambda2` |

The returned object is the **complete** prior: the residual-variance scale is
resolved, and for `mniw` the sum-of-coefficients / dummy-initial-observation
pseudo-observations are already folded into `Mu`/`V`/`Omega`/`DoF`. Nothing is
patched in later.

```matlab
prior = minnesotabvarm(3, 4, Y);                          % conjugate, defaults
prior = minnesotabvarm(3, 4, Y, lambda4=1, lambda5=1);    % + Sims-Zha dummies
prior = minnesotabvarm(3, 4, Y, Type="inw", lambda2=0.5); % asymmetric shrinkage
post  = estimate(prior);                                  % sample already stored
```

---

## 3. Hyperparameters

| Name | Meaning | Families | Default |
|---|---|---|---|
| `lambda1` | overall Minnesota tightness (paper λ) | all | `0.2` |
| `lambda2` | cross-variable relative tightness | `inw`, `normal` only | `1` (pinned in `mniw`) |
| `lambda3` | lag-decay exponent; variance ∝ 1/ℓ^(2·lambda3) | all | `1` |
| `lambda4` | sum-of-coefficients tightness (paper μ) | `mniw` only | `Inf` (off) |
| `lambda5` | dummy-initial-observation tightness (paper δ) | `mniw` only | `Inf` (off) |
| `Vc` | prior variance of constant / trend coefficients | all | `1e4` |
| `PriorMean` | prior mean on the own first lag | all | `1` |
| `Psi` | residual-variance scale (IW diagonal) | all | `"exact"` |

`lambda4 = lambda5 = Inf` means "dummy off" — the switch is `isfinite`, so an
`Inf` contributes no pseudo-observation rows (a zero row would still inflate the
degrees of freedom). The inverse-Wishart degrees of freedom are fixed at `m+2`
(the minimal proper choice with a defined prior mean); you do not set them.

`Psi` accepts a 1-by-`m` numeric vector, or the estimator names `"exact"`
(ARIMA-based, the historical default) / `"conditional"` (VARM-based, faster).
Inside a tuning loop, pass `Psi` numerically — the string form refits an AR per
series on every call.

---

## 4. Choosing hyperparameter values

- Litterman starting point: `PriorMean=1`, `lambda1≈0.2`, `lambda3≈1`, `Vc` large
  (`1e4`–`1e6`). Turn on `lambda4`/`lambda5` (e.g. `1`) for persistent macro data.
- Differenced / stationary series: set `PriorMean=0` for those series.
- Larger systems: tighter `lambda1` (more shrinkage) usually forecasts better.
- Rather than guess, tune by marginal likelihood — see §5.

---

## 5. Tuning by marginal likelihood (`glp` + `hyperprior`)

`glp` maximises the analytic marginal likelihood (plus any hyperpriors) over the
free hyperparameters, for the conjugate `mniw` family only. Each hyperparameter
is in one of three states, so the free/fixed split is implicit in what you pass:

| Value passed | State |
|---|---|
| scalar | **fixed** at that value |
| `[lo hi]` | **free**, flat prior on `[lo hi]` |
| `hyperprior(...)` | **free**, that log density added to the objective |

```matlab
psi = svar.estimateResidualVariances(Y, p, Method="conditional");
[prior, info] = glp(m, p, Y, ...
    lambda1 = hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5]), ...
    lambda4 = hyperprior("Gamma", 1,   1,   Bounds=[1e-4 50]), ...
    lambda5 = hyperprior("Gamma", 1,   1,   Bounds=[1e-4 50]), ...
    Vc      = 1e6, ...
    Psi     = psi);
```

`hyperprior(dist, mode, sd, Bounds=…)` takes an interpretable `(mode, sd)` pair
for `"Gamma"`, `"InverseGamma"`, or `"Beta"`. Some reference densities have no
mean or variance and cannot be reached from `(mode, sd)` — e.g. the
inverse-Gamma with shape = scale = `0.02^2` that GLP put on the residual scale.
Pass those natively:

```matlab
psiPrior = arrayfun(@(v) hyperprior("InverseGamma", 0.02^2, 0.02^2, ...
    Parameterization="native", X0=v, Bounds=[v/100 v*100]), psi);
[prior, info] = glp(m, p, Y, lambda1=..., Psi=psiPrior);   % free psi, one prior per series
```

`Psi` in `glp` accepts: a name (`"exact"`/`"conditional"`, resolved once), a
1-by-`m` vector (fixed), a 2-by-`m` `[lower; upper]` matrix (free, flat), a scalar
`hyperprior` (broadcast), or a 1-by-`m` `hyperprior` array.

**Integrating over the hyperparameters.** Ask for a third output plus `NumDraws`
to run GLP's Metropolis step instead of only maximising:

```matlab
[prior, info, chain] = glp(m, p, Y, lambda1=..., Psi=psiPrior, ...
    NumDraws=10000, BurnIn=10000, ProposalScale=0.8);
```

`chain` holds per-draw hyperparameters (`chain.lambda1`, `chain.Psi`, …) and the
VAR parameters they imply (`chain.Coefficients` `m·k × m × NumDraws`,
`chain.Sigma` `m × m × NumDraws`), plus `chain.AcceptanceRate`. Tune
`ProposalScale` for ~0.2–0.3 acceptance (`2.38/√d` for `d` free hyperparameters
is a good start). Set `rng` before the call. The chain draws feed stage 5
directly. `info.Hessian` / `info.ProposalCovariance` expose the proposal geometry.

---

## 6. The object-state rule (who carries the sample)

Some objects embed the estimation sample; some do not. Get this wrong and you
either double-count the data or condition on the wrong sample.

**Carries the sample** — call `estimate` / `simulate` / `forecast` with **no**
`Y`:
- `minnesotabvarm` outputs (the scale and dummies derive from the stored `Y`),
- the `glp` prior output,
- any posterior returned by `estimate`.

Passing a *different* `Y` to these errors on purpose.

**Data-free** — you **must** pass `Y`:
- `diffusebvarm`, `conjugatebvarm`, `semiconjugatebvarm`, `normalbvarm`,
- `weakbvarm`, `uniformirbvarm`.

```matlab
% data-free prior
[Coeff, Sigma] = simulate(weakbvarm(m,p), Y, NumDraws=Nd);

% sample-carrying posterior
post = estimate(minnesotabvarm(m,p,Y), Display="off");
[Coeff, Sigma] = simulate(post, NumDraws=Nd);          % no Y
```

For semiconjugate / `inw` workflows, `estimate` returns an `empiricalbvarm`
summary; when the next step needs raw coefficient/covariance draws, prefer
`simulate(prior, Y, NumDraws=Nd)`.

---

## 7. SSVS — what is and isn't available

- **Not** available for the `*bvarm` VAR objects (built-in or toolbox).
- **Available for Bayesian linear regression**: `mixsemiconjugateblm` performs
  SSVS; `customblm` allows arbitrary user-defined priors.
- For a VAR: either approximate with heavy `semiconjugatebvarm` / `inw`
  shrinkage, or estimate each equation with `mixsemiconjugateblm` and assemble
  the system. Be explicit about which; the BVAR objects do not do SSVS.
