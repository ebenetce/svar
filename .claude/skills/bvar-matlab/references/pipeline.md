# BVAR Pipeline Patterns Reference

Copy-paste MATLAB patterns for stages 2–5 that need more than one line. `Y` is a
T-by-m numeric matrix or timetable; `prior` / `Posterior` are BVAR objects from
stage 1–2. Follow the **object-state rule** (`references/priors.md` §6): pass `Y`
only to data-free priors.

## Table of contents
1. Lag-order selection
2. Forecast with predictive intervals and fan chart
3. Full predictive paths with `simsmooth`
4. IRF / FEVD credible bands (per-draw loop)
5. Hierarchical bands from a `glp` chain
6. FEVD plotting
7. Historical decomposition

---

## 1. Lag-order selection

The BVAR constructors do not pick lags. Use a frequentist `varm` pass for
AIC/BIC as a starting point, then build the BVAR at the chosen order:

```matlab
maxP = 8; aic = zeros(maxP,1); bic = zeros(maxP,1);
for p = 1:maxP
    EstMdl = estimate(varm(size(Y,2), p), Y);
    s = summarize(EstMdl);
    aic(p) = s.AIC; bic(p) = s.BIC;
end
[~,pAIC] = min(aic); [~,pBIC] = min(bic);
```

For BVARs a generous lag length (4 for quarterly, 12 for monthly) plus shrinkage
is common; the prior handles the over-parameterisation.

---

## 2. Forecast with predictive intervals and fan chart

`forecast` on the posterior object returns the mean path and its standard
deviations (data-free priors need `Y`; sample-carrying objects do not):

```matlab
fh = 12;
[YF, YFStd] = forecast(Posterior, fh, Y);      % fh-by-m each

z = norminv(0.95);                             % 90% band; 0.975 for 95%
lower = YF - z*YFStd;
upper = YF + z*YFStd;

tiledlayout(size(YF,2), 1);
for j = 1:size(YF,2)
    nexttile
    histT = (1:size(Y,1))';  fT = size(Y,1) + (1:fh)';
    plot(histT, Y(:,j), "k"); hold on
    fill([fT; flipud(fT)], [lower(:,j); flipud(upper(:,j))], ...
         [0.8 0.8 0.95], EdgeColor="none");
    plot(fT, YF(:,j), "b", LineWidth=1.2); hold off
    title(Posterior.SeriesNames(j))
end
```

---

## 3. Full predictive paths with `simsmooth`

For full predictive densities (not just mean±sd), append NaNs over the horizon
and let the simulation smoother fill them:

```matlab
numPaths = 1000; fh = 12;
YAug = [Y; NaN(fh, size(Y,2))];
[~, ~, NaNDraws] = simsmooth(Posterior, YAug, NumDraws=numPaths);
YSim = reshape(NaNDraws, fh, size(Y,2), numPaths);   % horizon x series x draw
qLo  = quantile(YSim, 0.05, 3);
qHi  = quantile(YSim, 0.95, 3);
```

Pinning some future values in `YAug` (leaving the rest NaN) gives a poor man's
**conditional forecast** until the dedicated API lands.

---

## 4. IRF / FEVD credible bands (per-draw loop)

The canonical structural loop: simulate reduced-form draws, convert each with
`svar.varmFromCoefficients`, identify, and apply. `svar.companionPower` is the
expensive impact-independent piece — compute it per draw and pass it to both
`svar.irf` and `svar.fevd`.

```matlab
H = 20; numDraws = 1000; shock = 1;   % shock index in SeriesNames order
rng(1);
[Coeff, Sigma] = simulate(Posterior, NumDraws=numDraws);
m = Posterior.NumSeries;

irfDraws  = zeros(H+1, m, numDraws);
fevdDraws = zeros(H+1, m, numDraws);
for d = 1:numDraws
    % simulate returns each draw as a vector; reshape to (m·k)-by-m for the converter
    Mdl    = svar.varmFromCoefficients(Posterior, reshape(Coeff(:,d), [], m), Sigma(:,:,d));
    impact = chol(Sigma(:,:,d), "lower");     % recursive ID in SeriesNames order
    Phi    = svar.companionPower(Mdl, H);
    ir     = svar.irf(Mdl,  impact, H, Phi=Phi);   % (H+1) x series x shock
    fe     = svar.fevd(Mdl, impact, H, Phi=Phi);
    irfDraws(:,:,d)  = ir(:,:,shock);
    fevdDraws(:,:,d) = fe(:,:,shock);
end

irfMed = median(irfDraws, 3);
irfLo  = quantile(irfDraws, 0.16, 3);   % 68% band
irfHi  = quantile(irfDraws, 0.84, 3);
```

`svar.irf` / `svar.fevd` already return `horizon × response × shock` — no
`permute`/`reshape`. To flip a shock's sign, negate that column of `impact`.

---

## 5. Hierarchical bands from a `glp` chain

When the prior came from `glp(..., NumDraws=n)`, the chain already carries VAR
parameter draws that integrate over the hyperparameters — skip `simulate`:

```matlab
H = 20; shock = 7;                     % e.g. monetary-policy shock, ordered last
numDraws = size(chain.Sigma, 3);
m = size(chain.Sigma, 1);
irfDraws = zeros(H+1, m, numDraws);
for d = 1:numDraws
    Mdl    = svar.varmFromCoefficients(prior, chain.Coefficients(:,:,d), chain.Sigma(:,:,d));
    impact = chol(chain.Sigma(:,:,d), "lower");   % recursive ID; funds rate last
    irfDraws(:,:,d) = svar.irf(Mdl, impact(:,shock), H);
end
irfBand = quantile(irfDraws, [0.16 0.84], 3);
```

Bands from the chain are wider than bands built at a single (modal)
hyperparameter value — that gap is the point of the hierarchical treatment.

---

## 6. FEVD plotting

```matlab
Decomp = median(fevdDraws, 3);   % (H+1) x series, share of `shock` (from §4)

% Or the full variable x shock decomposition at the posterior mode:
Mdl    = bvar2var(Posterior);
impact = chol(Mdl.Covariance, "lower");
fe     = svar.fevd(Mdl, impact, H);      % (H+1) x series x shock

tiledlayout(m, 1);
for i = 1:m
    nexttile
    area(squeeze(fe(:,i,:)));            % shares of each shock for variable i
    title("FEVD of " + Posterior.SeriesNames(i));
    if i == 1, legend(Posterior.SeriesNames, Location="eastoutside"); end
end
```

---

## 7. Historical decomposition

`historicalDecomposition` takes a **BVAR object or a `varm`** (a BVAR is
converted internally via `bvar2var`), a **square** `impact`, and the sample:

```matlab
impact = chol(bvar2var(Posterior).Covariance, "lower");   % recursive ID; comment the ordering
HD     = historicalDecomposition(Posterior, impact, Y);

% HD.Contributions : (T-p) x numVariables x numShocks (cumulative shock contributions)
% HD.Total         : (T-p) x numVariables (sum over selected shocks)
% HD.StructuralShocks, HD.VariableNames, HD.ShockNames, ...
```

Restrict the output with `VariableIndices` / `ShockIndices`, supply presample /
predictors with `Y0=` / `X=`, or pass precomputed `Residuals=`. Plot
`HD.Contributions(:,i,:)` as a stacked bar/area per variable `i` against the
demeaned series.
