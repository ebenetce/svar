%DEMOFIVESTEPPIPELINE Five-step BVAR/SVAR demo for Banbura et al. (2008).
%
%   This is the intentionally simple version of the replication. It uses the
%   official JAE data prepared by PREPAREBANBURA2008DATA and shows the same
%   modeling surface for the SMALL, CEE and MEDIUM systems:
%
%       Step 1: prior
%       Step 2: estimate
%       Step 3: simulate
%       Step 4: identify
%       Step 5: apply
%
%   The core five calls stay in the main script. Helper functions only make
%   the figures and print compact FEVD summaries.
%
%   See also MINNESOTABVARM, ESTIMATE, SIMULATE, BVAR2VAR, SVAR.IRF,
%   SVAR.FEVD.

clear
clc

thisFolder = fileparts(mfilename("fullpath"));
addpath(thisFolder);

load(fullfile(thisFolder, "banbura2008Data.mat"), "tt", "priorMean");
specs = banbura2008ModelSpecifications();

modelNames = ["SMALL" "CEE" "MEDIUM"];
targetVariables = ["EMPL" "CPI" "FFR"];
numLags = 13;
horizon = 48;
numDraws = 300;
rng(0, "twister");

% The paper's structural exercise estimates 1961-2002.
sample = tt(timerange(datetime(1961,1,1), datetime(2002,12,1), "closed"), :);

fprintf("Banbura, Giannone and Reichlin (2008): five-step BVAR/SVAR demo\n");
fprintf("Sample: %s to %s. Lags: %d. Posterior draws: %d.\n\n", ...
    string(sample.Time(1), "yyyy-MM"), string(sample.Time(end), "yyyy-MM"), ...
    numLags, numDraws);

for modelName = modelNames
    spec = specs.(modelName);
    Y = sample{:, spec.Variables};
    thisPriorMean = localSelectPriorMean(tt, priorMean, spec.Variables);
    numSeries = size(Y, 2);

    fprintf("=== %s (%d variables) ===\n", modelName, numSeries);

    %% Step 1: prior
    % This is the paper's conjugate Minnesota prior with the sum-of-
    % coefficients dummy switched on: tau = 10*lambda.
    priorMdl = minnesotabvarm(numSeries, numLags, Y, ...
        lambda1 = spec.lambda1, ...
        lambda3 = 1, ...
        lambda4 = spec.lambda4, ...
        Psi = "conditional", ...
        PriorMean = thisPriorMean, ...
        SeriesNames = spec.Variables);

    %% Step 2: estimate
    % Estimate is inherited from Econometrics Toolbox's conjugatebvarm
    % surface; the repo's Minnesota class stores Y, so no sample argument is
    % needed here.
    posterior = estimate(priorMdl, Display = "off");

    %% Step 3: simulate
    % Draw reduced-form coefficients and covariance matrices from the
    % posterior. These are used for credible bands.
    [coeffDraws, sigmaDraws] = simulate(posterior, NumDraws = numDraws);

    %% Step 4: identify
    % Recursive monetary-policy shock: variables are ordered slow-moving
    % first, then FFR, then fast-moving variables. We scale the FFR shock to
    % a 100 bp impact response.
    pointVAR = bvar2var(posterior);
    impact = chol(pointVAR.Covariance, "lower");
    policyShock = impact(:, spec.PolicyIndex) ...
        ./ impact(spec.PolicyIndex, spec.PolicyIndex);

    %% Step 5: apply
    % Apply the structural shock with this repo's array-returning helpers.
    pointIRF = svar.irf(pointVAR, policyShock, horizon);
    pointFEVD = svar.fevd(pointVAR, impact, horizon);

    irfDraws = zeros(horizon + 1, numSeries, numDraws);
    for draw = 1:numDraws
        drawVAR = svar.varmFromCoefficients(posterior, ...
            reshape(coeffDraws(:, draw), [], numSeries), sigmaDraws(:,:,draw));
        drawImpact = chol(sigmaDraws(:,:,draw), "lower");
        drawPolicyShock = drawImpact(:, spec.PolicyIndex) ...
            ./ drawImpact(spec.PolicyIndex, spec.PolicyIndex);
        irfDraws(:,:,draw) = svar.irf(drawVAR, drawPolicyShock, horizon);
    end

    localPrintFEVD(modelName, pointFEVD, spec, targetVariables, [1 3 6 12 24 36 48]);
    localPlotIRF(modelName, pointIRF, irfDraws, spec.Variables, targetVariables);
end

function values = localSelectPriorMean(tt, priorMean, names)
allNames = string(tt.Properties.VariableNames);
idx = arrayfun(@(name) find(allNames == name, 1), names);
values = priorMean(idx);
end

function localPrintFEVD(~, fevd, spec, targetVariables, horizons)
fprintf("FEVD share of the monetary-policy shock, percent\n");
for variable = targetVariables
    variableIndex = find(spec.Variables == variable, 1);
    share = 100*fevd(horizons + 1, variableIndex, spec.PolicyIndex);
    fprintf("  %-5s %s\n", variable, mat2str(round(share(:)', 1)));
end
fprintf("\n");
end

function localPlotIRF(modelName, pointIRF, irfDraws, names, targetVariables)
horizons = 0:(size(pointIRF, 1) - 1);
targetIndex = arrayfun(@(name) find(names == name, 1), targetVariables);
bands = quantile(irfDraws(:, targetIndex, :), [0.05 0.16 0.84 0.95], 3);

figure(Color = "w", Name = modelName + " five-step demo");
layout = tiledlayout(numel(targetVariables), 1, ...
    TileSpacing = "compact", Padding = "compact");

for j = 1:numel(targetVariables)
    ax = nexttile(layout);
    fill(ax, [horizons fliplr(horizons)], ...
        [bands(:,j,1); flipud(bands(:,j,4))]', [0.88 0.90 0.95], ...
        EdgeColor = "none");
    hold(ax, "on");
    fill(ax, [horizons fliplr(horizons)], ...
        [bands(:,j,2); flipud(bands(:,j,3))]', [0.72 0.79 0.92], ...
        EdgeColor = "none");
    plot(ax, horizons, pointIRF(:, targetIndex(j)), "k", LineWidth = 1.2);
    yline(ax, 0, ":");
    xlim(ax, [0 horizons(end)]);
    title(ax, targetVariables(j));
    box(ax, "on");
end

title(layout, modelName + ": response to a 100 bp federal funds rate shock");
end
