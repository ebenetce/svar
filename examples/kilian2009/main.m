%MAIN Reproduce Kilian (2009), "Not All Oil Price Shocks Are Alike".
%
%   American Economic Review 99(3), 1053-1069. A three-variable structural
%   VAR of the global crude-oil market, identified recursively.
%
%   The point of this file is that a faithful reproduction is short. Kilian's
%   prior is a diffuse Normal-Inverse-Wishart, so the model IS one line -
%   diffusebvarm - and the rest is the standard BVAR pipeline: estimate,
%   simulate, identify (Cholesky), apply (irf / fevd / historical
%   decomposition). Everything runs on the Econometrics Toolbox object surface.
%
%   System, ordered most- to least-exogenous (Kilian's recursive ordering):
%     1  global crude oil production growth    -> oil SUPPLY shock
%     2  index of global real economic activity -> aggregate DEMAND shock
%     3  real price of oil                      -> oil-specific DEMAND shock
%
%   Data is fetched to TEMPDIR at run time, so no data file is committed.
%   The vintage differs from Kilian's 2009 archive (his real-activity index
%   has since been revised), so the numbers are close in spirit, not identical.
%
%   See also DIFFUSEBVARM, SVAR.IRF, SVAR.FEVD, HISTORICALDECOMPOSITION.

clear
clc

numLags    = 24;      % Kilian: 24 monthly lags
irfHorizon = 15;      % months
numDraws   = 500;     % posterior draws for the credible bands
rng(0, "twister")

%% Data
data = localFetchKilianData();
Y    = [100*[NaN; diff(data.Oil_Prod)], data.Kilian_Index, 100*data.RAC];
keep = data.Data <= datetime(2007,12,1);        % Kilian's 1973M1-2007M12 window
Y    = Y(keep, :);
Y    = Y(2:end, :);                             % drop the leading diff NaN
dates = data.Data(keep);
dates = dates(numLags+2:end);                   % dates aligned to the residuals

seriesNames = ["Oil production growth", "Real activity", "Real oil price"];
shockNames  = ["Oil supply", "Aggregate demand", "Oil-specific demand"];

%% Prior and estimate — Kilian's diffuse NIW prior
posterior = estimate(diffusebvarm(3, numLags, SeriesNames=seriesNames), ...
    Y, Display="off");

%% Identify recursively and apply at the posterior mean
% Lower-triangular impact in the ordering above; the supply column is negated
% so a positive shock is a supply DISRUPTION (production falls, price rises).
pointVAR = bvar2var(posterior);
impact   = localRecursiveImpact(pointVAR.Covariance);

irf  = svar.irf(pointVAR,  impact, irfHorizon);   % (H+1) x variable x shock
fevd = svar.fevd(pointVAR, impact, 12);           % share at h = 12
HD   = historicalDecomposition(posterior, impact, Y, VariableIndices=3);

%% Credible bands: identify and apply per posterior draw
[coeffDraws, sigmaDraws] = simulate(posterior, NumDraws=numDraws);
irfDraws = zeros(irfHorizon+1, 3, 3, numDraws);
for d = 1:numDraws
    drawVAR = svar.varmFromCoefficients(posterior, ...
        reshape(coeffDraws(:,d), [], 3), sigmaDraws(:,:,d));
    irfDraws(:,:,:,d) = svar.irf(drawVAR, ...
        localRecursiveImpact(sigmaDraws(:,:,d)), irfHorizon);
end
band = quantile(irfDraws, [0.16 0.84], 4);        % 68% band

%% Impulse responses (Kilian, Figure 3)
horizons = 0:irfHorizon;
figure(Color="w");
layout = tiledlayout(3, 3, TileSpacing="compact", Padding="compact");
for variable = 1:3
    for shock = 1:3
        ax = nexttile(layout);
        fill(ax, [horizons, fliplr(horizons)], ...
            [band(:,variable,shock,1); flipud(band(:,variable,shock,2))]', ...
            [0.85 0.88 0.95], EdgeColor="none");
        hold(ax, "on");
        plot(ax, horizons, irf(:,variable,shock), "k", LineWidth=1.1);
        yline(ax, 0, ":");
        xlim(ax, [0 irfHorizon]);
        if variable == 1, title(ax, shockNames(shock)); end
        if shock == 1, ylabel(ax, seriesNames(variable)); end
    end
end
title(layout, "Kilian (2009): responses to one-s.d. structural shocks");

%% Variance decomposition of the real oil price at h = 12
fprintf("Share of the real-oil-price forecast-error variance at h = 12:\n");
for shock = 1:3
    fprintf("  %-20s %5.1f%%\n", shockNames(shock), 100*fevd(end, 3, shock));
end

%% Historical decomposition of the real oil price
figure(Color="w");
area(dates, squeeze(HD.Contributions(:, 1, :)));
hold on
plot(dates, Y(numLags+1:end, 3) - mean(Y(numLags+1:end, 3)), "k", LineWidth=1);
legend([shockNames, "demeaned real oil price"], Location="northwest");
title("Historical decomposition of the real oil price");
xlabel("date"); ylabel("percent"); box on

%% Local functions
function impact = localRecursiveImpact(Sigma)
%LOCALRECURSIVEIMPACT Cholesky impact with Kilian's supply-shock sign flip.
impact = chol(Sigma, "lower");
impact(:, 1) = -impact(:, 1);
end

function data = localFetchKilianData()
%LOCALFETCHKILIANDATA Load the oil-market series, caching in TEMPDIR.
%   Oil_Prod (log global oil production), Kilian_Index (global real activity),
%   RAC (log real refiner acquisition cost = the real price of oil).
cacheFile = fullfile(tempdir, "kilian2009_Data_oil_1.xlsx");
if ~isfile(cacheFile)
    base = "https://raw.githubusercontent.com/mauep2025/Global-Oil-Market/";
    fetched = false;
    for branch = ["master", "main"]
        try
            websave(cacheFile, base + branch + "/Data_oil_1.xlsx");
            fetched = true;
            break
        catch
        end
    end
    if ~fetched
        error("kilian2009:downloadFailed", ...
            "Could not download the data. Place Data_oil_1.xlsx at %s.", cacheFile);
    end
end
data      = readtable(cacheFile);
data.Data = datetime(data.Data);
end
