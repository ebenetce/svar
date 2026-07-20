# BVAR Pipeline Patterns Reference

Copy-paste MATLAB patterns for pipeline steps that need more than one line. All
examples assume `Y` is a T-by-m numeric matrix or timetable and `PriorMdl` is a
BVAR prior object that can be conditioned on `Y`.

## Table of contents
1. Lag-order selection
2. Forecast with predictive intervals and fan chart
3. Constructing `varm` objects from BVAR draws
4. IRF with posterior credible bands
5. FEVD extraction and plotting
6. Historical decomposition caveat

---

## 1. Lag-order selection

The BVAR constructors do not pick lags. Use a frequentist `varm` pass for
AIC/BIC as a starting point, then build the BVAR at the chosen order:

```matlab
maxP = 8; aic = zeros(maxP,1); bic = zeros(maxP,1);
for p = 1:maxP
    Mdl = varm(size(Y,2), p);
    EstMdl = estimate(Mdl, Y);
    results = summarize(EstMdl);
    aic(p) = results.AIC; bic(p) = results.BIC;
end
[~,pAIC] = min(aic); [~,pBIC] = min(bic);
```

For BVARs, a generous lag length such as 4 for quarterly data or 12 for monthly
data plus shrinkage is common; shrinkage handles much of the
over-parameterization.

---

## 2. Forecast with predictive intervals and fan chart

`forecast(PriorMdl, fh, Y)` is valid: MATLAB conditions the BVAR prior on the
supplied sample.

```matlab
fh = 12;
[YF, YFStd] = forecast(PriorMdl, fh, Y);       % fh-by-m

z = norminv(0.95);                 % 90% band; use 0.975 for 95%
lower = YF - z*YFStd;
upper = YF + z*YFStd;

m = size(YF,2);
figure;
for j = 1:m
    subplot(m,1,j);
    hold on;
    histT = (1:size(Y,1))';
    fT    = size(Y,1) + (1:fh)';
    plot(histT, Y(:,j), 'k');
    fill([fT; flipud(fT)], [lower(:,j); flipud(upper(:,j))], ...
         [0.8 0.8 0.95], EdgeColor="none");
    plot(fT, YF(:,j), 'b', LineWidth=1.2);
    title(PriorMdl.SeriesNames(j)); hold off;
end
```

For full predictive paths, append NaNs for the forecast horizon and use
`simsmooth`; reshape the NaN draws by horizon and variable:

```matlab
numPaths = 1000;
YAug = [Y; NaN(fh, size(Y,2))];
[~,~,NaNDraws,YMean,YStd] = simsmooth(PriorMdl, YAug, NumDraws=numPaths);
YSim = reshape(NaNDraws, fh, size(Y,2), numPaths);
qLo = quantile(YSim, 0.05, 3);
qHi = quantile(YSim, 0.95, 3);
```

---

## 3. Constructing `varm` objects from BVAR draws

`irf` and `fevd` are `varm` methods. Do not call them directly on `*bvarm`
objects. For BVAR structural analysis, draw coefficients and innovation
covariances from the data-conditioned BVAR, convert each draw to a `varm`, and
then call `irf` / `fevd`.

The coefficient ordering in `CoeffDraws` follows the object's coefficient map,
not a universal hand-written ordering. Always inspect:

```matlab
S = summarize(PriorMdl);
disp(S.CoeffMap)
```

For standard BVARs with constants and AR coefficients, use `CoeffMap` to unpack
draws:

```matlab
function Mdl = bvarDrawToVARM(PriorMdl, coeff, Sigma)
    S = summarize(PriorMdl);
    map = string(S.CoeffMap);
    m = PriorMdl.NumSeries;
    p = PriorMdl.P;

    Mdl = varm(m, p);
    Mdl.SeriesNames = PriorMdl.SeriesNames;
    Mdl.Covariance = Sigma;

    for lag = 1:p
        A = zeros(m);
        for response = 1:m
            for predictor = 1:m
                key = "AR{" + lag + "}(" + response + "," + predictor + ")";
                idx = find(map == key, 1);
                if ~isempty(idx)
                    A(response,predictor) = coeff(idx);
                end
            end
        end
        Mdl.AR{lag} = A;
    end

    c = zeros(m,1);
    hasConstant = false;
    for response = 1:m
        idx = find(map == "Constant(" + response + ")", 1);
        if ~isempty(idx)
            c(response) = coeff(idx);
            hasConstant = true;
        end
    end
    if hasConstant
        Mdl.Constant = c;
    end
end
```

Extend the helper for trend and exogenous `Beta` terms when those are enabled.

---

## 4. IRF with posterior credible bands

For credible bands, draw from the data-conditioned BVAR and recompute the VAR
IRF per draw. `irf` returns `H`-by-`m`-by-`m` arrays, horizon by response by
shock, orthogonalized by Cholesky in `SeriesNames` order by default.

```matlab
H = 20; numDraws = 1000;
[Coeff, Sigma] = simulate(PriorMdl, Y, NumDraws=numDraws);
m = PriorMdl.NumSeries;
irfDraws = zeros(H, m, m, numDraws);
for d = 1:numDraws
    Mdl_d = bvarDrawToVARM(PriorMdl, Coeff(:,d), Sigma(:,:,d));
    irfDraws(:,:,:,d) = irf(Mdl_d, NumObs=H);
end
irfMed = median(irfDraws, 4);
irfLo  = quantile(irfDraws, 0.16, 4);  % 68% band
irfHi  = quantile(irfDraws, 0.84, 4);
```

Sign convention: to flip a shock's sign, multiply that shock's IRF slice by -1
because VARs are linear.

---

## 5. FEVD extraction and plotting

```matlab
H = 20; numDraws = 1000;
[Coeff, Sigma] = simulate(PriorMdl, Y, NumDraws=numDraws);
m = PriorMdl.NumSeries;
fevdDraws = zeros(H, m, m, numDraws);
for d = 1:numDraws
    Mdl_d = bvarDrawToVARM(PriorMdl, Coeff(:,d), Sigma(:,:,d));
    fevdDraws(:,:,:,d) = fevd(Mdl_d, NumObs=H);
end
Decomp = median(fevdDraws, 4);  % H-by-m-by-m: horizon x variable x shock

figure;
for i = 1:size(Decomp,2)
    subplot(size(Decomp,2),1,i);
    area(squeeze(Decomp(:,i,:)));
    title("FEVD of " + PriorMdl.SeriesNames(i));
    if i==1, legend(PriorMdl.SeriesNames, Location="eastoutside"); end
end
```

---

## 6. Historical decomposition caveat

There is no single `histdecomp` call for BVAR objects in the Econometrics
Toolbox. To build one, convert each BVAR draw to a `varm`, infer reduced-form
residuals, orthogonalize them with the chosen structural impact matrix such as
the Cholesky factor of `Sigma`, and convolve past structural shocks with the
draw-specific IRFs. Tell the user this is a manual construction and confirm they
want it before writing the loop.
