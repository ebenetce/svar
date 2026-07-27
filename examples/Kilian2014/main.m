%MAIN Replicate the reduced-form Kilian-Murphy (2014) oil-market VAR.
%
%   Kilian and Murphy, "The Role of Inventories and Speculative Trading in
%   the Global Market for Crude Oil", Journal of Applied Econometrics 29(3),
%   454-478.
%
%   This script reproduces the parts of the paper that do not require the
%   paper's sign/zero/admissibility restrictions: the data set, sample,
%   diffuse Gaussian-inverse Wishart reduced-form posterior, two years of
%   monthly lags, seasonal dummies, posterior simulation, and reduced-form
%   diagnostics.
%
%   The paper's named shocks are set-identified by sign restrictions,
%   elasticity bounds, dynamic restrictions, and Haar rotations. Those steps
%   are intentionally omitted here. The IRFs, FEVDs, and historical
%   decompositions below use a recursive Cholesky impact matrix as a clearly
%   labelled diagnostic, not as the paper's benchmark structural model.
%
%   See also DIFFUSEBVARM, SVAR.IRF, SVAR.FEVD, HISTORICALDECOMPOSITION.

clear
clc

thisFolder = fileparts(mfilename("fullpath"));
repoRoot   = fileparts(fileparts(thisFolder));
proj       = openProject(fullfile(repoRoot, "SVAR.prj"));

fprintf("Project: %s\n", proj.Name);

%% Settings
numLags       = 24;     % two years of monthly lags
irfHorizon    = 20;     % Figure 1 reports responses through 20 months
fevdHorizon   = 60;
numDraws      = 112;    % paper reports 112 reduced-form posterior draws
randomSeed    = 0;
priceVariable = 3;

seriesNames = [ ...
    "Oil production growth", ...
    "Real activity", ...
    "Real price of oil", ...
    "Inventory change"];

recursiveShockNames = [ ...
    "Production innovation", ...
    "Activity innovation", ...
    "Oil-price innovation", ...
    "Inventory innovation"];

%% Data
data = localLoadKilianMurphyData(thisFolder);
Y     = data.Y;
dates = data.Dates;

% The sample starts in February, so the default reference category
% (period 12) is the following January.
X = svar.seasonalDummies(size(Y, 1), 12);

fprintf("Sample: %s to %s, %d monthly observations\n", ...
    string(dates(1), "yyyy-MM"), string(dates(end), "yyyy-MM"), size(Y, 1));
fprintf("Variables: %s\n\n", strjoin(seriesNames, ", "));

%% Reduced-form posterior
% Kilian and Murphy use a diffuse Gaussian-inverse Wishart prior for the
% reduced-form parameters. Seasonal dummies enter as exogenous predictors.
prior = diffusebvarm(numel(seriesNames), numLags, ...
    SeriesNames   = seriesNames, ...
    NumPredictors = size(X, 2));

posterior = estimate(prior, Y, X = X, Display = "off");
pointVAR  = bvar2var(posterior);
residuals = infer(pointVAR, Y, X = X);

fprintf("Reduced-form posterior: %s, %d lags, %d seasonal dummies\n", ...
    class(posterior), posterior.P, posterior.NumPredictors);
fprintf("Effective residual sample: %d observations\n\n", size(residuals, 1));

disp("Posterior mean innovation covariance:");
disp(array2table(pointVAR.Covariance, ...
    VariableNames = seriesNames, RowNames = seriesNames));

%% Recursive diagnostic impact matrix
% Recursive ordering: oil production growth, real activity, real oil price,
% inventory change. This is only a diagnostic normalization; it is not the
% sign-restricted Kilian-Murphy identification.
impact = localRecursiveImpact(pointVAR.Covariance, priceVariable);

disp("Recursive impact matrix, price-normalized where possible:");
disp(array2table(impact, VariableNames = recursiveShockNames, ...
    RowNames = seriesNames));

fprintf("Recursive first-shock production/price impact ratio: %.3f\n\n", ...
    impact(1, 1)/impact(3, 1));

%% Figure 1-style impulse responses
irf = svar.irf(pointVAR, impact, irfHorizon);
irf = localAccumulateFlowVariables(irf);

rng(randomSeed, "twister");
[coeffDraws, sigmaDraws] = simulate(posterior, NumDraws = numDraws);

irfDraws = zeros(irfHorizon + 1, numel(seriesNames), ...
    numel(seriesNames), numDraws);
for draw = 1:numDraws
    drawVAR = svar.varmFromCoefficients(posterior, ...
        reshape(coeffDraws(:, draw), [], numel(seriesNames)), ...
        sigmaDraws(:, :, draw));
    drawImpact = localRecursiveImpact(sigmaDraws(:, :, draw), priceVariable);
    drawIRF    = svar.irf(drawVAR, drawImpact, irfHorizon);
    irfDraws(:, :, :, draw) = localAccumulateFlowVariables(drawIRF);
end
irfBand = quantile(irfDraws, [0.16 0.84], 4);

localPlotIRFs(irf, irfBand, seriesNames, recursiveShockNames(1:3), ...
    "Kilian-Murphy (2014): recursive Figure 1-style diagnostic");

%% FEVD diagnostic
fevd = svar.fevd(pointVAR, impact, fevdHorizon);
localPrintFEVD(fevd, fevdHorizon, seriesNames, recursiveShockNames);

%% Historical decompositions, recursive diagnostic
residualDates = dates((numLags + 1):end);
HD = historicalDecomposition(pointVAR, impact, Y, ...
    X = X, VariableIndices = [3 4], ShockIndices = 1:3);

eventDates = [ ...
    datetime(1978, 9, 1), ...
    datetime(1980, 9, 1), ...
    datetime(1985,12, 1), ...
    datetime(1990, 8, 1), ...
    datetime(1997, 7, 1), ...
    datetime(2002,11, 1)];

localPlotLongHD(residualDates, HD.Contributions(:, 1, :), ...
    recursiveShockNames(1:3), eventDates);

episodes = [ ...
    struct(Name = "Persian Gulf War, 1990/91", ...
        Start = datetime(1990, 7, 1), Stop = datetime(1991, 2, 1))
    struct(Name = "Iranian Revolution, 1978/79", ...
        Start = datetime(1978,10, 1), Stop = datetime(1980, 1, 1))
    struct(Name = "Iran-Iraq War outbreak, 1980", ...
        Start = datetime(1980, 9, 1), Stop = datetime(1981, 6, 1))
    struct(Name = "Collapse of OPEC, 1986", ...
        Start = datetime(1986, 1, 1), Stop = datetime(1986, 6, 1))
    struct(Name = "Venezuelan crisis and Iraq War, 2002/03", ...
        Start = datetime(2002,11, 1), Stop = datetime(2003, 5, 1))];

for episode = episodes.'
    localPlotEpisodeHD(residualDates, HD.Contributions, ...
        recursiveShockNames(1:3), episode);
end

%% Replication gaps
fprintf("\nReplication gaps intentionally left out:\n");
fprintf("  * Table 1 impact sign restrictions and the dynamic sign restrictions.\n");
fprintf("  * Supply and demand impact-elasticity bounds used to select admissible models.\n");
fprintf("  * Haar rotations: the paper draws 5 million rotations and keeps 14 admissible models.\n");
fprintf("  * Paper Figure 1 bands, Figures 2-7 structural decompositions, and Table 2 elasticities.\n");
fprintf("    The plots above are recursive diagnostics using the same reduced-form posterior.\n");

%% Local functions
function data = localLoadKilianMurphyData(thisFolder)
load(fullfile(thisFolder, "kmData.mat"), "kmData");
load(fullfile(thisFolder, "worldprod.mat"), "worldprod");

dates = (datetime(1973, 2, 1):calmonths(1):datetime(2009, 8, 1)).';
if size(kmData, 1) ~= numel(dates)
    error("Kilian2014:DataSizeMismatch", ...
        "kmData must have one row for each month from 1973M2 to 2009M8.");
end
if numel(worldprod) ~= size(kmData, 1) + 1
    error("Kilian2014:WorldProductionSizeMismatch", ...
        "worldprod should run from 1973M1 to 2009M8.");
end

data = struct();
data.Y = kmData;
data.Dates = dates;
data.WorldProduction = worldprod;
end

function impact = localRecursiveImpact(Sigma, priceVariable)
impact = chol(Sigma, "lower");

% Match the paper's normalization where possible: reported shocks move the
% real price of oil up on impact. The last Cholesky shock has zero impact on
% the price by construction and is left with its positive diagonal sign.
for shock = 1:min(priceVariable, size(impact, 2))
    if impact(priceVariable, shock) < 0
        impact(:, shock) = -impact(:, shock);
    end
end
end

function irf = localAccumulateFlowVariables(irf)
% Oil production is a growth rate and inventories are changes.
irf(:, [1 4], :) = cumsum(irf(:, [1 4], :), 1);
end

function localPlotIRFs(irf, band, seriesNames, shockNames, plotTitle)
horizons = 0:(size(irf, 1) - 1);
numVariables = numel(seriesNames);
numShocks = numel(shockNames);

figure(Color = "w");
layout = tiledlayout(numVariables, numShocks, ...
    TileSpacing = "compact", Padding = "compact");

for variable = 1:numVariables
    for shock = 1:numShocks
        ax = nexttile(layout);
        fill(ax, [horizons, fliplr(horizons)], ...
            [band(:, variable, shock, 1); ...
             flipud(band(:, variable, shock, 2))].', ...
            [0.86 0.89 0.95], EdgeColor = "none");
        hold(ax, "on");
        plot(ax, horizons, irf(:, variable, shock), ...
            "k", LineWidth = 1.1);
        yline(ax, 0, ":", Color = [0.45 0.45 0.45]);
        xlim(ax, [horizons(1) horizons(end)]);
        box(ax, "on");

        if variable == 1
            title(ax, shockNames(shock));
        end
        if shock == 1
            ylabel(ax, seriesNames(variable));
        end
        if variable == numVariables
            xlabel(ax, "Months");
        end
    end
end

title(layout, plotTitle);
end

function localPrintFEVD(fevd, horizon, seriesNames, shockNames)
for h = [12 horizon]
    fprintf("Recursive FEVD shares at h = %d months:\n", h);
    rows = 100*squeeze(fevd(h + 1, [3 4], :));
    disp(array2table(rows, VariableNames = shockNames, ...
        RowNames = seriesNames([3 4])));
end
end

function localPlotLongHD(dates, priceContributions, shockNames, eventDates)
plotStart = datetime(1978, 6, 1);
keep = dates >= plotStart;
hdates = dates(keep);
priceContributions = squeeze(priceContributions(keep, 1, :));

figure(Color = "w");
layout = tiledlayout(numel(shockNames), 1, ...
    TileSpacing = "compact", Padding = "compact");

for shock = 1:numel(shockNames)
    ax = nexttile(layout);
    plot(ax, hdates, priceContributions(:, shock), "k", LineWidth = 1.1);
    hold(ax, "on");
    yline(ax, 0, ":", Color = [0.45 0.45 0.45]);
    for eventDate = eventDates
        xline(ax, eventDate, ":", Color = [0.6 0.6 0.6]);
    end
    ylabel(ax, shockNames(shock));
    box(ax, "on");
end

xlabel(layout, "Date");
title(layout, "Recursive Figure 2-style real-oil-price decomposition");
end

function localPlotEpisodeHD(dates, contributions, shockNames, episode)
keep = dates >= episode.Start & dates <= episode.Stop;
if ~any(keep)
    return
end

episodeDates = dates(keep);
episodeContributions = contributions(keep, :, :);

figure(Color = "w");
layout = tiledlayout(2, 1, TileSpacing = "compact", Padding = "compact");

for variable = 1:2
    ax = nexttile(layout);
    plot(ax, episodeDates, squeeze(episodeContributions(:, variable, :)), ...
        LineWidth = 1.1);
    yline(ax, 0, ":", Color = [0.45 0.45 0.45]);
    box(ax, "on");
    if variable == 1
        ylabel(ax, "Real price");
        legend(ax, shockNames, Location = "best");
    else
        ylabel(ax, "Inventories");
    end
end

xlabel(layout, "Date");
title(layout, "Recursive event decomposition: " + episode.Name);
end
