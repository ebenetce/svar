# BVAR Priors, Hyperparameters & Extension Reference

Primary path is the **direct constructors** (`conjugatebvarm`,
`semiconjugatebvarm`, `normalbvarm`, `diffusebvarm`), which take the full prior
specification. `bayesvarm` is only the Minnesota-shortcut wrapper — avoid unless
the user asks for it.

## Table of contents
1. The coefficient vector λ and the Mu / V layout
2. Building Minnesota Mu and V by hand
3. Constructor cheat-sheet and hyperparameters
4. Choosing hyperparameter values
5. Subclassing the prior objects (extension recipes)
6. Sims–Zha dummy-observation prior
7. SSVS — what is and isn't available
8. Construct-as-posterior idiom (reduced-form SVAR) + manual-sampler Rosetta stone

---

## 1. The coefficient vector λ and the Mu / V layout

For an m-variable VAR(p) with constant c, optional trend δ, and r exogenous
predictors Β, the per-equation coefficient block is
`[Φ1 Φ2 … Φp c δ Β]'`, an `(m*p + 1 + hasTrend + r)`-by-`m` matrix Λ.
The estimable vector is `λ = vec(Λ')` (stack equations). Let
`k = m*p + IncludeConstant + IncludeTrend + NumPredictors` be the number of
coefficients per equation; then `λ` has `m*k` elements.

- `Mu`  — prior mean of `λ`, an `(m*k)`-by-1 vector.
- `V`   — prior covariance **factor** of the coefficients (for conjugate models
  the full conditional cov is Σ⊗V; for semiconjugate, V is the coefficient cov).
- `Omega` — inverse-Wishart scale matrix Ω (m-by-m, PD).
- `DoF`   — inverse-Wishart degrees of freedom ν (proper if ν>m−1; finite mean
  if ν>m+1).
- `Sigma` — fixed innovations covariance (only `normalbvarm`).

Inspect the exact ordering on any constructed object via `PriorMdl.Mu`,
`PriorMdl.V`, and the derived `PriorMdl.AR`, `.Constant`, `.Trend`, `.Beta`.

---

## 2. Building Minnesota Mu and V by hand

This reproduces what `bayesvarm` does internally, but under your control.

Minnesota rules:
- Prior mean: first own lag = `center` (1 for levels / random-walk prior, 0 for
  stationary/differenced data); all other coefficients mean 0. Fill `Mu`
  accordingly.
- Prior variance of the coefficient on lag ℓ of variable j in equation i:
  - own lag (i=j):  `(lambda / ℓ^decay)^2`
  - cross lag (i≠j): `(lambda * cross / ℓ^decay)^2 * (σ_i^2 / σ_j^2)`
  - exogenous/deterministic: a large value `varianceX` (near-flat).
  where `σ_i^2` are residual variances from univariate AR fits (the scale
  factors), `lambda` is overall tightness, `cross∈(0,1]` downweights other
  variables, `decay≥0` is the lag-decay rate.

Sketch:
```matlab
function Mu = buildMinnesotaMu(m,p,center,hasConst,hasTrend,r)
    k  = m*p + hasConst + hasTrend + r;
    Mu = zeros(m*k,1);
    for i = 1:m                      % equation i
        base = (i-1)*k;
        Mu(base + i) = center;       % coeff on own first lag
    end
end

function V = buildMinnesotaV(m,p,lambda,cross,decay,sig,varX,hasConst,hasTrend,r)
    k = m*p + hasConst + hasTrend + r;
    v = zeros(m*k,1);
    for i = 1:m
        base = (i-1)*k; col = 1;
        for ell = 1:p
            for j = 1:m
                if i==j, s = (lambda/ell^decay)^2;
                else,    s = (lambda*cross/ell^decay)^2*(sig(i)/sig(j)); end
                v(base+col) = s; col = col+1;
            end
        end
        for d = 1:(hasConst+hasTrend+r)   % deterministic / exogenous
            v(base+col) = varX; col = col+1;
        end
    end
    V = diag(v);
end
```
(Indexing assumes the per-equation block orders lags then deterministics; verify
against `PriorMdl.Mu` ordering on a small example before trusting it.)

---

## 3. Constructor cheat-sheet

| Constructor | Prior | estimate → | Set |
|---|---|---|---|
| `diffusebvarm(m,p)` | flat | analytic `conjugatebvarm` | nothing |
| `conjugatebvarm(m,p,...)` | Matrix-Normal-IW (Σ⊗V) | analytic `conjugatebvarm` | Mu,V,Omega,DoF |
| `semiconjugatebvarm(m,p,...)` | independent Normal-IW | Gibbs → `empiricalbvarm` | Mu,V,Omega,DoF |
| `normalbvarm(m,p,...)` | normal coeffs, Σ fixed | analytic | Mu,V,Sigma |

Conjugacy constraint for `conjugatebvarm`: the coefficient covariance must be
proportional across equations (Σ⊗V Kronecker structure → symmetric own/cross
treatment). If you need asymmetric own-vs-cross variances, you cannot use the
plain conjugate object — use `semiconjugatebvarm` or subclass (see §5).

---

## 4. Choosing hyperparameter values

- Litterman starting point: `center=1`, `lambda≈0.2`, `cross≈0.5`, `decay≈1–2`,
  `varianceX` large (≈1e6), `DoF=m+2`, `Omega=diag(σ²)`.
- Differenced/stationary data: `center=0`.
- Larger systems: tighter `lambda` (more shrinkage) usually forecasts better.
- Tune `lambda` (optionally `cross`) by grid/axial search maximizing marginal
  likelihood (closed-form for the conjugate prior) or minimizing
  pseudo-out-of-sample loss. Loop over constructor calls; the toolbox doesn't
  automate it.

---

## 5. Subclassing the prior objects (extension recipes)

The `*bvarm` objects are classes; subclass the closest one and override in the
constructor (and optionally `estimate`/`forecast`/`simulate`).

- **Asymmetric conjugate prior** (own ≠ cross variances, but keep analytics):
  subclass `conjugatebvarm`; override `V` to the asymmetric structure. Note this
  breaks strict Σ⊗V conjugacy — document that the analytic posterior is then an
  approximation, or route `estimate` through a Gibbs step.
- **Sum-of-coefficients / Sims–Zha dummy prior baked in**: subclass
  `conjugatebvarm`; override `estimate` to append the dummy observations (see §6)
  to the data, then call the superclass `estimate`.
- **Hierarchical hyperprior on λ (e.g. random tightness)**: subclass
  `semiconjugatebvarm`; override `estimate` to add a Gibbs block that resamples
  the hyperparameter, updating `V` each sweep before delegating.
- **Fully custom prior density**: the regression-family analogue is `customblm`
  (user-declared log-prior function). There is no `custombvarm`, so for a VAR
  encode the custom prior in a subclass with an overridden sampler, or build it
  equation-by-equation with `customblm`.

Skeleton:
```matlab
classdef sumOfCoeffsBVARM < conjugatebvarm
    properties
        Mu_sc (1,1) double = 1   % sum-of-coefficients tightness
    end
    methods
        function obj = sumOfCoeffsBVARM(m,p,varargin)
            obj = obj@conjugatebvarm(m,p,varargin{:});
        end
        function Post = estimate(obj, Y, varargin)
            [Yd, ~] = buildSimsZhaDummies(obj, Y);   % §6
            Post = estimate@conjugatebvarm(obj, [Y; Yd], varargin{:});
        end
    end
end
```

---

## 6. Sims–Zha dummy-observation prior

Implement by data augmentation, then estimate a conjugate model on the stacked
data. Blocks:
- **Sum-of-coefficients** (persistence): rows built from presample means of each
  series, scaled by tightness μ; pushes toward unit-root behavior.
- **Co-persistence / dummy-initial-observation**: a single row tying the
  constant to the common stochastic trend, scaled by δ.

Procedure: compute presample series means → construct `Yd`,`Xd` per Sims–Zha
(1998) / Bańbura, Giannone & Reichlin (2010) with tightness (λ, μ, δ) → stack
`[Y;Yd]`,`[X;Xd]` → estimate conjugate. Confirm (λ, μ, δ) with the user and
write the dummy construction explicitly; don't approximate it with plain
Minnesota args and call it Sims–Zha.

---

## 7. SSVS — what is and isn't available

- **Not** available for the `*bvarm` VAR objects.
- **Available for Bayesian linear regression**: `mixsemiconjugateblm` performs
  SSVS; `customblm` allows arbitrary user-defined priors.
- For a VAR: either approximate with heavy `semiconjugatebvarm` shrinkage, or
  estimate each equation with `mixsemiconjugateblm` and assemble the system.
  Be explicit about which; the BVAR objects do not do SSVS.

---

## 8. Reduced-form SVAR draws: object-state rule + diffuse OLS shortcut

Reduced-form parameter draws for external identification (sign restrictions,
Haar rotations) follow the same API rule for all built-in BVAR priors
(`diffusebvarm`, `conjugatebvarm`, `semiconjugatebvarm`, `normalbvarm`) and for
subclasses:

**(a) Simulate from a prior object with data supplied** (valid and usually the
right route when the object is a pure prior that still needs to be conditioned
on `Y`, especially for semiconjugate models):
```matlab
[CoeffDraws, SigmaDraws] = simulate(PriorMdl, Y, NumDraws=Nd);
```

**(b) Estimate or construct a posterior, then simulate the posterior object**
(valid when the posterior object supports `simulate` and already carries the
estimation sample or posterior moments):
```matlab
PosteriorMdl = estimate(PriorMdl, Y, Display="off");
[CoeffDraws, SigmaDraws] = simulate(PosteriorMdl, NumDraws=Nd);
```

Do **not** pass `Y` to `simulate` when the object already represents a posterior
that embeds the estimation sample. That includes manually constructed posterior
objects and custom/subclassed objects that embed `Y`, dummy observations,
presample data, or transformed estimation data. Passing a new `Y` asks MATLAB to
condition the supplied object on that data, which can double-count the sample or
silently condition on a different sample than the one embedded in the object.
For semiconjugate workflows, remember that `estimate` returns an `empiricalbvarm`
summary object; if the next step needs raw coefficient and covariance draws,
prefer `simulate(PriorMdl, Y, NumDraws=Nd)` or `simsmooth(PriorMdl, Y, ...)`.

### Diffuse-prior OLS NIW construct-as-posterior shortcut

The OLS recipe below is **only** for the diffuse/flat-prior posterior, or for
intentionally matching legacy hand-rolled SVAR code that uses the same NIW
moments. Do not use it for Minnesota, semiconjugate, normal, or custom priors
unless their posterior algebra has been explicitly derived to match these
moments.

Under a diffuse prior the conjugate posterior is closed-form, so you can fit the
VAR by OLS, compute the posterior NIW moments yourself, load them straight into
a `conjugatebvarm` object, and simulate that object. This is useful when you
want a specific degrees-of-freedom correction or transparent replication-grade
moments.

```matlab
% OLS fit (VARX with seasonal dummies as predictors)
Mdl0  = varm(K, p);
Est   = estimate(Mdl0, y, X=Predictors);
SIGMA = Est.Covariance * (t-p) / (t-p - K*p - 1);   % df-corrected (NOT the MLE Sigma)
betaOLS = [ [Est.AR{:}], Est.Constant, Est.Beta ];  % see ordering caveat below

% Diffuse-prior posterior NIW moments
T     = t - p;
Vfac  = inv(X*X');        % coefficient covariance factor (Sigma ⊗ V)
Omega = T * SIGMA;        % inverse-Wishart scale
DoF   = T;                % inverse-Wishart degrees of freedom
Mu    = betaOLS(:);       % posterior mean = OLS coefficients

PostMdl = conjugatebvarm(K, p, IncludeConstant=true, NumPredictors=size(Predictors,2), ...
                         Mu=Mu, V=Vfac, Omega=Omega, DoF=DoF);
rng(1);
[CoeffDraws, SigmaDraws] = simulate(PostMdl, NumDraws=n1);  % draw from constructed posterior
% CoeffDraws: (numcoeff)-by-n1 ;  SigmaDraws: K-by-K-by-n1
```

For a standard diffuse reduced-form posterior, `diffusebvarm` + `estimate` +
`simulate(PosteriorMdl, NumDraws=Nd)` is valid. Use the OLS NIW construction
when you need the exact hand-derived moments, the chosen df correction, or a
one-to-one replacement for a legacy NIW sampler.

### Rosetta stone: hand-rolled NIW sampler ↔ `conjugatebvarm` properties

Legacy SVAR code often hand-codes the NIW posterior draw. The object reproduces
it exactly. Mapping the manual loop to object properties:

| Manual sampler step | Object equivalent |
|---|---|
| `RANTR = randn(T,q)/sqrt(T)*chol(inv(SIGMA)); SIGMAr = inv(RANTR'*RANTR)` — Wishart/Bartlett draw, then invert ⇒ Σ⁽ʳ⁾ ~ IW | `Omega = T*SIGMA`, `DoF = T` (the IW component) |
| `vecAr = Bvec + chol(kron(pXX, SIGMAr))' * randn(...)` — conditional-normal coeff draw, mean = OLS coeffs, cov = Σ⁽ʳ⁾ ⊗ (X'X)⁻¹ | `Mu = OLS coeffs`, `V = inv(X*X')` (the conditional-normal component) |
| Draw Σ first, then coefficients given Σ | exactly the NIW factorization `simulate` performs |

So the entire hand-rolled loop collapses to the `conjugatebvarm` construction
above plus `simulate(PostMdl, NumDraws=Nd)`. Recognizing this pattern in
someone's legacy code and offering the two-line toolbox equivalent is a
high-value move.

### Coefficient-ordering caveat (do not assume)

The order of coefficients in `Mu` / `betaOLS` is **construction-dependent** and
must match how `X` (and thus `V = inv(X*X')`) is laid out:
- Some code stacks **AR lags first, deterministics last** (`Mu = [AR(:); det(:)]`).
- Other code stacks **deterministics first, then AR lags** (e.g. `X = [const; dummies; lagged Y]`).
These are not interchangeable: `Mu`, `V`, and `X` must all use the same ordering,
or the draws are silently wrong. **Always verify** against `PriorMdl.Mu` and the
actual `X` layout on a small example before trusting the construction — print
`Est.AR{:}`, `Est.Constant`, `Est.Beta` and confirm the stacking by hand.
