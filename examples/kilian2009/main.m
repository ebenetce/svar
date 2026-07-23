%MAIN Reproduce Kilian (2009) with the SVAR toolbox — recursive identification.
%
%   "Not All Oil Price Shocks Are Alike: Disentangling Demand and Supply
%   Shocks in the Crude Oil Market", American Economic Review 99(3), 1053-1069.
%
%   Kilian identifies three structural shocks in the global crude-oil market
%   with a RECURSIVE (Cholesky, lower-triangular) scheme — no sign or zero
%   restrictions — so it is a clean exercise of exactly what this toolbox does
%   today: the Cholesky identify stage, and all three structural apply
%   functions, SVAR.IRF, SVAR.FEVD, and HISTORICALDECOMPOSITION.
%
%   The five-stage pipeline
%   -----------------------
%     1 prior     weakbvarm (improper NIW; posterior mean = OLS)
%     2 estimate  estimate(prior, Y)
%     3 simulate  simulate(post, NumDraws) for credible bands
%     4 identify  impact = chol(Sigma, "lower")  in Kilian's ordering
%     5 apply     svar.irf / svar.fevd / historicalDecomposition
%
%   The 3-variable system, ordered most-to-least exogenous:
%     1  global crude oil production growth   (oil SUPPLY shock)
%     2  index of global real economic activity (aggregate DEMAND shock)
%     3  real price of oil                    (oil-specific DEMAND shock)
%
%   Data & reproducibility
%   ----------------------
%   The three series are fetched at run time from a public replication archive
%   and cached in TEMPDIR — nothing is written into the repository, so the only
%   file this example commits is main.m. The vintage differs from Kilian's 2009
%   archive (his global real-activity index has since been revised), so the
%   point estimates are close in spirit, not identical to the published
%   figures. The PRECISE correctness check is at the end: on this exact sample
%   the weak-prior posterior mean equals OLS, and svar.irf equals the built-in
%   orthogonalized irf — both to machine precision, independent of vintage.
%
%   See also WEAKBVARM, SVAR.VARMFROMCOEFFICIENTS, SVAR.IRF, SVAR.FEVD,
%   HISTORICALDECOMPOSITION.

clear
clc

%% Settings
numLags     = 24;                 % Kilian: 24 monthly lags
irfHorizon  = 15;                 % months
numDraws    = 500;                % posterior draws for credible bands
randomSeed  = 0;
band        = [0.16 0.84];        % 68% credible band
sampleEnd   = datetime(2007,12,1);% Kilian's original window ends 2007M12

%% Data — fetched at run time, cached in TEMPDIR (never committed)
dataTable = localFetchKilianData();

oilProdGrowth = 100*[NaN; diff(dataTable.Oil_Prod)];  % % change in oil production
realActivity  = dataTable.Kilian_Index;               % global real activity index
realOilPrice  = 100*dataTable.RAC;                     % log real oil price (x100)

keep = dataTable.Data <= sampleEnd;
Y    = [oilProdGrowth, realActivity, realOilPrice];
Y    = Y(keep, :);
Y    = Y(2:end, :);                                   % drop the leading diff NaN
dates = dataTable.Data(keep);
dates = dates(2:end);

seriesNames = ["OilProdGrowth", "RealActivity", "RealOilPrice"];
shockNames  = ["Oil supply", "Aggregate demand", "Oil-specific demand"];
numSeries   = size(Y, 2);

fprintf("Kilian (2009) system: %d obs x %d vars, %s to %s\n\n", ...
    size(Y,1), numSeries, string(dates(1),"MMM-yyyy"), string(dates(end),"MMM-yyyy"));

%% Stages 1-2 — prior and estimate
% weakbvarm is an improper reduced-form NIW prior; its posterior mean equals
% the OLS estimate, so this reproduces Kilian's frequentist VAR while flowing
% through the Bayesian object surface.
prior     = weakbvarm(numSeries, numLags, SeriesNames = seriesNames);
posterior = estimate(prior, Y, Display = "off");

modeCoefficients = reshape(posterior.Mu, [], numSeries);
modeSigma        = posterior.Omega/(posterior.DoF + numSeries + 1);
modeVAR          = svar.varmFromCoefficients(posterior, modeCoefficients, modeSigma);

%% Stage 4 — identify (recursive / Cholesky, Kilian's ordering)
% Lower-triangular impact in SeriesNames order: oil production reacts to none
% within the month, real activity to supply only, the oil price to both. The
% oil-supply column is negated so a positive shock is a supply DISRUPTION
% (production falls, price rises) — Kilian's sign convention.
modeImpact         = chol(modeSigma, "lower");
modeImpact(:, 1)   = -modeImpact(:, 1);

%% Stage 3 — draws for credible bands
rng(randomSeed, "twister");
[coefficientDraws, covarianceDraws] = simulate(posterior, NumDraws = numDraws);

irfDraws = zeros(irfHorizon + 1, numSeries, numSeries, numDraws);
for draw = 1:numDraws
    sigma   = covarianceDraws(:,:,draw);
    impact  = chol(sigma, "lower");
    impact(:, 1) = -impact(:, 1);
    drawVAR = svar.varmFromCoefficients(posterior, ...
        reshape(coefficientDraws(:,draw), [], numSeries), sigma);
    irfDraws(:,:,:,draw) = svar.irf(drawVAR, impact, irfHorizon);
end

%% Stage 5 — structural impulse responses (Kilian, Figure 3)
modeIRF = svar.irf(modeVAR, modeImpact, irfHorizon);      % (H+1) x variable x shock
irfBand = quantile(irfDraws, band, 4);
horizons = 0:irfHorizon;

figure(Color = "w");
layout = tiledlayout(numSeries, numSeries, TileSpacing = "compact", Padding = "compact");
for variable = 1:numSeries
    for shock = 1:numSeries
        ax = nexttile(layout);
        fill(ax, [horizons, fliplr(horizons)], ...
            [irfBand(:,variable,shock,1); flipud(irfBand(:,variable,shock,2))]', ...
            [0.85 0.88 0.95], EdgeColor = "none");
        hold(ax, "on");
        plot(ax, horizons, modeIRF(:,variable,shock), "k", LineWidth = 1.1);
        yline(ax, 0, ":", Color = [0.5 0.5 0.5]);
        xlim(ax, [0 irfHorizon]);
        if variable == 1, title(ax, shockNames(shock)); end
        if shock == 1, ylabel(ax, seriesNames(variable)); end
    end
end
title(layout, "Kilian (2009): responses to one-s.d. structural shocks");

%% Stage 5 — FEVD of the real price of oil
fevd12 = svar.fevd(modeVAR, modeImpact, 12);      % (13) x variable x shock
fprintf("Share of the real-oil-price forecast-error variance at h = 12:\n");
for shock = 1:numSeries
    fprintf("  %-20s %5.1f%%\n", shockNames(shock), 100*fevd12(end, 3, shock));
end

%% Stage 5 — historical decomposition of the real price of oil
HD = historicalDecomposition(posterior, modeImpact, Y, VariableIndices = 3);

figure(Color = "w");
area(dates(numLags+1:end), squeeze(HD.Contributions(:, 1, :)));
hold on
plot(dates(numLags+1:end), Y(numLags+1:end, 3) - mean(Y(numLags+1:end, 3)), ...
    "k", LineWidth = 1);
title("Historical decomposition of the real oil price");
legend([shockNames, "demeaned real oil price"], Location = "northwest");
xlabel("date"); ylabel("percent"); box on

%% Correctness cross-checks (exact, independent of the data vintage)
olsVAR = estimate(varm(numSeries, numLags), Y);
maxCoeffGap = max(abs(cell2mat(modeVAR.AR) - cell2mat(olsVAR.AR)), [], "all");

builtinIRF = permute(irf(modeVAR, NumObs = irfHorizon + 1, ...
    Method = "orthogonalized"), [1 3 2]);
plainImpact = chol(modeSigma, "lower");
maxIRFGap = max(abs(svar.irf(modeVAR, plainImpact, irfHorizon) - builtinIRF), [], "all");

fprintf("\nCorrectness checks:\n");
fprintf("  max |weak posterior-mean AR - OLS AR|        = %.2e\n", maxCoeffGap);
fprintf("  max |svar.irf - builtin orthogonalized irf|  = %.2e\n", maxIRFGap);

%% Local functions
function dataTable = localFetchKilianData()
%LOCALFETCHKILIANDATA Load the oil-market series, caching in TEMPDIR.
%   Fetched from a public replication repository so this example needs no
%   committed data file. Columns used: Oil_Prod (log global oil production),
%   Kilian_Index (global real activity), RAC (log real refiner acquisition
%   cost = the real price of oil).
cacheFile = fullfile(tempdir, "kilian2009_Data_oil_1.xlsx");

if ~isfile(cacheFile)
    base = "https://raw.githubusercontent.com/mauep2025/Global-Oil-Market/";
    sources = base + ["master"; "main"] + "/Data_oil_1.xlsx";
    fetched = false;
    for source = sources'
        try
            websave(cacheFile, source);
            fetched = true;
            break
        catch
        end
    end
    if ~fetched
        error("kilian2009:downloadFailed", ...
            "Could not download the oil-market data. Check the network, or " + ...
            "place Data_oil_1.xlsx at %s manually.", cacheFile);
    end
end

dataTable      = readtable(cacheFile);
dataTable.Data = datetime(dataTable.Data);
end
