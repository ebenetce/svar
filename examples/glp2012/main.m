%main Native driver for the Giannone-Lenza-Primiceri (2012) example.

clear
clc;

glpNumDraws = 20000;
numLags = 5;
forecastHorizons = 1:8;
irfHorizon = 20;
policyShock = 7;

load('data.mat', 'tt')

numSeries = width(tt);

% Select Prior
spec = minnesotaSpec( "mniw",...
    lambda1 = hyperprior('Gamma', 0.2, 0.4, Bounds = [1e-4, 5]), ...
    lambda3 = 1, ...
    lambda4 = hyperprior('Gamma', 1, 1, Bounds = [1e-4, 50]), ...
    lambda5 = hyperprior('Gamma', 1, 1, Bounds = [1e-4, 50]), ...
    Vc = 10e6);

Psi = estimateResidualVariances(tt{:,:}, 1, Method = "conditional");
Psi = arrayfun(@(x) hyperprior("InverseGamma", 0.02^2, 0.02^2, ...
    X0 = x, Bounds = [1/100, 100]*x), Psi);

priorMdl = glp(numSeries, numLags, tt, Psi, ...
    Spec = spec);

% Estimate Posterior
posterior = estimate(prior, Y);

%% 
glpOptions = struct();
glpOptions.Spec = spec;
glpOptions.PsiScale = [1/100, 100];
glpOptions.PsiPriorMode = 0.02^2;
glpOptions.PsiPriorSD = 0.02^2;
glpOptions.PsiLags = 1;
glpOptions.PsiMethod = "conditional";
glpOptions.IncludeConstant = true;
glpOptions.IncludeTrend = false;
glpOptions.SeriesNames = seriesNames;
glpOptions.OptimOptions = optimoptions("fmincon", ...
    Display = "final", ...
    FiniteDifferenceStepSize = 1e-4, ...
    FunctionTolerance = 1e-12, ...
    StepTolerance = 1e-12, ...
    ConstraintTolerance = 1e-12, ...
    MaxIterations = 1000, ...
    MaxFunctionEvaluations = 5000);

defaultRun = localRunGLPModel(responses, numLags, forecastHorizons, ...
    irfHorizon, policyShock, glpNumDraws, glpRandomSeed, glpOptions, ...
    "GLP all-priors");

localPrintHyperparameters(defaultRun);
localPrintForecastTable(defaultRun, time(end), "Default GLP prior");

if glpMakePlots
    densityFigure = localPlotPredictiveDensity(defaultRun, time(end));
    if glpExportPlots
        localSaveFigure(densityFigure, ...
            fullfile(resultsFolder, "default_predictive_density.png"));
    end

    defaultIRFFigure = localPlotIRFs(defaultRun, displayNames, ...
        "Default GLP prior");
    if glpExportPlots
        localSaveFigure(defaultIRFFigure, ...
            fullfile(resultsFolder, "default_monetary_policy_irfs.png"));
    end
end

economicChecks = localEconomicChecks(defaultRun);
localPrintEconomicChecks(economicChecks);

if glpRunLegacyAudit
    legacyAudit = localLegacyAudit(responses, numLags, forecastHorizons, ...
        defaultRun, thisFolder);
    localPrintLegacyAudit(legacyAudit);
else
    legacyAudit = [];
end

if glpRunMinnesotaOnly
    mnOnlyOptions = glpOptions;
    mnOnlyOptions.Spec.lambda4 = Inf;
    mnOnlyOptions.Spec.lambda5 = Inf;

    mnOnlyRun = localRunGLPModel(responses, numLags, forecastHorizons, ...
        irfHorizon, policyShock, glpNumDraws, glpRandomSeed + 1, ...
        mnOnlyOptions, "Minnesota-only");

    fprintf("\n");
    localPrintHyperparameters(mnOnlyRun);
    localPrintForecastTable(mnOnlyRun, time(end), "Minnesota-only prior");

    if glpMakePlots
        mnOnlyIRFFigure = localPlotIRFs(mnOnlyRun, displayNames, ...
            "Minnesota-only prior");
        if glpExportPlots
            localSaveFigure(mnOnlyIRFFigure, ...
                fullfile(resultsFolder, "minnesota_only_monetary_policy_irfs.png"));
        end
    end
else
    mnOnlyRun = [];
end

results = struct( ...
    "Project", proj.Name, ...
    "DataFile", dataFile, ...
    "Series", displayNames, ...
    "Time", time, ...
    "NumLags", numLags, ...
    "ForecastHorizons", forecastHorizons, ...
    "IRFHorizon", irfHorizon, ...
    "PolicyShock", policyShock, ...
    "Default", defaultRun, ...
    "MinnesotaOnly", mnOnlyRun, ...
    "EconomicChecks", economicChecks, ...
    "LegacyAudit", legacyAudit);

save(fullfile(resultsFolder, "glp2012_results.mat"), "results");
fprintf("\nSaved results to %s\n", ...
    fullfile(resultsFolder, "glp2012_results.mat"));

function coef = localGammaCoef(modeValue, sdValue)
    ratio = modeValue^2/sdValue^2;
    coef.k = (2 + ratio + sqrt((4 + ratio)*ratio))/2;
    coef.theta = sqrt(sdValue^2/coef.k);
end

function run = localRunGLPModel(Y, numLags, horizons, irfHorizon, ...
        policyShock, numDraws, randomSeed, options, label)
    numSeries = size(Y, 2);
    psi0 = estimateResidualVariances(Y, options.PsiLags, ...
        Method = options.PsiMethod);
    Psi = arrayfun(@(x) hyperprior("InverseGamma", ...
        options.PsiPriorMode, options.PsiPriorSD, ...
        X0 = x, Bounds = options.PsiScale*x), psi0);

    [prior, info] = glp(numSeries, numLags, Y, Psi, ...
        Spec = options.Spec, ...
        IncludeConstant = options.IncludeConstant, ...
        IncludeTrend = options.IncludeTrend, ...
        SeriesNames = options.SeriesNames, ...
        OptimOptions = options.OptimOptions);

    posterior = estimate(prior, Y, Display="off");
    modeCoefficients = reshape(posterior.Mu, [], numSeries);
    modeSigma = posterior.Omega/(posterior.DoF + numSeries + 1);

    modeForecast = localDeterministicForecast(Y, modeCoefficients, horizons);
    modeGDPGrowth = localAverageGDPGrowth(modeForecast(:, 1), ...
        Y(end, 1), horizons);

    modeVAR = varmFromCoefficients(posterior, modeCoefficients, modeSigma);
    modeImpact = chol(modeSigma, "lower");
    modeIRF = svar.irf(modeVAR, modeImpact(:, policyShock), irfHorizon - 1);

    rng(randomSeed, "twister");
    [coefficientDraws, covarianceDraws] = simulate(posterior, ...
        NumDraws=numDraws);

    forecastDraws = zeros(numel(horizons), numSeries, numDraws);
    irfDraws = zeros(irfHorizon, numSeries, numDraws);

    for draw = 1:numDraws
        coefficients = reshape(coefficientDraws(:, draw), [], numSeries);
        sigma = covarianceDraws(:, :, draw);

        forecastDraws(:, :, draw) = localStochasticForecast(Y, ...
            coefficients, sigma, horizons);

        drawVAR = varmFromCoefficients(posterior, coefficients, sigma);
        impact = chol(sigma, "lower");
        irfDraws(:, :, draw) = svar.irf(drawVAR, ...
            impact(:, policyShock), irfHorizon - 1);
    end

    gdpForecastDraws = squeeze(forecastDraws(:, 1, :));
    gdpGrowthDraws = localAverageGDPGrowth(gdpForecastDraws, ...
        Y(end, 1), horizons);

    sortedIRFDraws = sort(irfDraws, 3);
    bandIndex = max(1, min(numDraws, round([0.16, 0.50, 0.84]*numDraws)));
    irfBands = sortedIRFDraws(:, :, bandIndex);

    run = struct( ...
        "Label", label, ...
        "Prior", prior, ...
        "Posterior", posterior, ...
        "Info", info, ...
        "ModeCoefficients", modeCoefficients, ...
        "ModeSigma", modeSigma, ...
        "ModeForecast", modeForecast, ...
        "ModeGDPGrowth", modeGDPGrowth, ...
        "ForecastDraws", forecastDraws, ...
        "GDPGrowthDraws", gdpGrowthDraws, ...
        "ModeIRF", modeIRF, ...
        "IRFDraws", irfDraws, ...
        "IRFBands", irfBands, ...
        "NumDraws", numDraws, ...
        "RandomSeed", randomSeed);
end

function forecast = localDeterministicForecast(Y, coefficients, horizons)
    numSeries = size(Y, 2);
    numLags = (size(coefficients, 1) - 1)/numSeries;
    yPath = [Y; zeros(max(horizons), numSeries)];
    start = size(Y, 1);

    for step = 1:max(horizons)
        regressors = localLagRegressors(yPath, start + step, ...
            numLags, numSeries);
        yPath(start + step, :) = regressors*coefficients;
    end

    forecast = yPath(start + horizons, :);
end

function forecast = localStochasticForecast(Y, coefficients, sigma, horizons)
    numSeries = size(Y, 2);
    numLags = (size(coefficients, 1) - 1)/numSeries;
    yPath = [Y; zeros(max(horizons), numSeries)];
    shockChol = chol(sigma, "lower");
    start = size(Y, 1);

    for step = 1:max(horizons)
        regressors = localLagRegressors(yPath, start + step, ...
            numLags, numSeries);
        shock = shockChol*randn(numSeries, 1);
        yPath(start + step, :) = regressors*coefficients + shock.';
    end

    forecast = yPath(start + horizons, :);
end

function regressors = localLagRegressors(yPath, row, numLags, numSeries)
    regressors = zeros(1, numLags*numSeries + 1);
    for lag = 1:numLags
        cols = ((lag - 1)*numSeries + 1):(lag*numSeries);
        regressors(cols) = yPath(row - lag, :);
    end
    regressors(end) = 1;
end

function growth = localAverageGDPGrowth(gdpForecast, lastGDP, horizons)
    growth = (gdpForecast - lastGDP) .* (100 ./ horizons(:));
end

function localPrintHyperparameters(run)
    spec = run.Info.FinalSpec;
    fprintf("%s hyperparameters:\n", run.Label);
    fprintf("  lambda1 = %.6g\n", spec.lambda1);
    fprintf("  lambda3 = %.6g\n", spec.lambda3);
    fprintf("  lambda4 = %.6g\n", spec.lambda4);
    fprintf("  lambda5 = %.6g\n", spec.lambda5);
    fprintf("  psi     = %s\n", mat2str(run.Info.FinalPsi, 6));
    fprintf("  objective = %.6f, exitflag = %.0f\n", ...
        run.Info.Objective, run.Info.ExitFlag);
end

function localPrintForecastTable(run, forecastDate, label)
    horizonsToPrint = [1, 4, 8];
    drawStats = [median(run.GDPGrowthDraws(horizonsToPrint, :), 2), ...
        std(run.GDPGrowthDraws(horizonsToPrint, :), 0, 2)];
    tableData = [horizonsToPrint(:), run.ModeGDPGrowth(horizonsToPrint), ...
        drawStats];

    fprintf("\n%s\n", label);
    fprintf("Forecasting in %s (GDP annualized average growth)\n", ...
        localQuarterLabel(forecastDate, 2));
    fprintf(" Qrts ahead    Point       Median        Std\n");
    fprintf("            (at Mode)   predictive  predictive\n");
    fprintf(" %6.0f      %8.4f    %8.4f    %8.4f\n", tableData.');
end

function fig = localPlotPredictiveDensity(run, forecastDate)
    fig = figure(Color="w");
    chart = tiledlayout(2, 1, TileSpacing="compact", Padding="compact");
    horizons = [1, 4];

    for idx = 1:numel(horizons)
        horizon = horizons(idx);
        ax = nexttile(chart);
        histogram(ax, run.GDPGrowthDraws(horizon, :), 50, ...
            FaceColor=[0.25 0.45 0.68], EdgeColor="none");
        xline(ax, run.ModeGDPGrowth(horizon), Color=[0.1 0.1 0.1], ...
            LineWidth=1.25);
        title(ax, sprintf("Forecasting in %s: Avg GDP growth %d qrts ahead", ...
            localQuarterLabel(forecastDate, 2), horizon));
        xlabel(ax, "Annualized average growth");
        ylabel(ax, "Draw count");
        xlim(ax, [-15, 10]);
        box(ax, "on");
    end
end

function fig = localPlotIRFs(run, displayNames, label)
    horizons = 0:(size(run.ModeIRF, 1) - 1);
    fig = figure(Color="w");
    chart = tiledlayout(2, 4, TileSpacing="compact", Padding="compact");

    for variable = 1:numel(displayNames)
        ax = nexttile(chart);
        modeResponse = localScaleIRF(run.ModeIRF(:, variable), variable);
        bandResponse = localScaleIRF( ...
            squeeze(run.IRFBands(:, variable, :)), variable);

        plot(ax, horizons, modeResponse, Color=[0 0 0], LineWidth=1.1);
        hold(ax, "on");
        plot(ax, horizons, bandResponse(:, 1), "-.", Color=[0.75 0.1 0.1]);
        plot(ax, horizons, bandResponse(:, 2), "-.", Color=[0.75 0.1 0.1]);
        plot(ax, horizons, bandResponse(:, 3), "-.", Color=[0.75 0.1 0.1]);
        yline(ax, 0, Color=[0.45 0.45 0.45], LineStyle=":");
        xlim(ax, [horizons(1), horizons(end)]);
        title(ax, displayNames(variable));
        box(ax, "on");

        if variable == numel(displayNames)
            legend(ax, "IRF at mode", "16th quantile", ...
                "50th quantile", "84th quantile", ...
                Location="northeastoutside");
        end
    end

    title(chart, label + ": monetary-policy shock");
end

function checks = localEconomicChecks(run)
    scaledIRF = zeros(size(run.ModeIRF));
    for variable = 1:size(run.ModeIRF, 2)
        scaledIRF(:, variable) = localScaleIRF(run.ModeIRF(:, variable), ...
            variable);
    end

    checks = struct();
    checks.PolicyRateImpact = scaledIRF(1, 7);
    checks.AverageGDPResponseH1To8 = mean(scaledIRF(2:9, 1));
    checks.AveragePriceResponseH4To12 = mean(scaledIRF(5:13, 2));
    checks.GDPGrowthProbabilityBelowZero = mean( ...
        run.GDPGrowthDraws([1, 4, 8], :) < 0, 2);
end

function localPrintEconomicChecks(checks)
    fprintf("\nEconomic-behavior checks for default GLP prior:\n");
    fprintf("  Fed funds impact response: %.4f basis points\n", ...
        checks.PolicyRateImpact);
    fprintf("  Avg GDP response, horizons 1-8: %.4f percent\n", ...
        checks.AverageGDPResponseH1To8);
    fprintf("  Avg GDP deflator response, horizons 4-12: %.4f percent\n", ...
        checks.AveragePriceResponseH4To12);
    fprintf("  Pr(avg GDP growth < 0), h = 1/4/8: %s\n", ...
        mat2str(checks.GDPGrowthProbabilityBelowZero.', 4));
end

function audit = localLegacyAudit(Y, numLags, horizons, currentRun, thisFolder)
    legacyFolder = fullfile(thisFolder, "lenzaPrimiceri");
    oldPath = path;
    cleanup = onCleanup(@() path(oldPath));
    addpath(legacyFolder);

    legacyRun = localQuietLegacyGLP(Y, numLags);

    currentCoeffLegacyOrder = [currentRun.ModeCoefficients(end, :); ...
        currentRun.ModeCoefficients(1:end-1, :)];
    currentGDPGrowth = currentRun.ModeGDPGrowth;
    legacyGDPGrowth = localAverageGDPGrowth(legacyRun.postmax.forecast(:, 1), ...
        Y(end, 1), horizons);

    audit = struct();
    audit.Note = "Legacy GLP MAP benchmark, used as a health check rather than an exact target.";
    audit.RelativeCoefficientDifference = localRelativeNorm( ...
        currentCoeffLegacyOrder, legacyRun.postmax.betahat);
    audit.RelativeSigmaDifference = localRelativeNorm( ...
        currentRun.ModeSigma, legacyRun.postmax.sigmahat);
    audit.RelativeForecastDifference = localRelativeNorm( ...
        currentRun.ModeForecast, legacyRun.postmax.forecast);
    audit.CurrentHyperparameters = struct( ...
        "lambda1", currentRun.Info.FinalSpec.lambda1, ...
        "lambda4SumOfCoefficients", currentRun.Info.FinalSpec.lambda4, ...
        "lambda5DummyInitialObservation", currentRun.Info.FinalSpec.lambda5, ...
        "alphaEquivalent", 2*currentRun.Info.FinalSpec.lambda3, ...
        "logPosteriorObjective", -currentRun.Info.Objective);
    audit.LegacyHyperparameters = struct( ...
        "lambda", legacyRun.postmax.lambda, ...
        "miuSumOfCoefficients", legacyRun.postmax.miu, ...
        "thetaDummyInitialObservation", legacyRun.postmax.theta, ...
        "alpha", legacyRun.postmax.alpha, ...
        "logPost", legacyRun.postmax.logPost);
    audit.GDPGrowthHorizons = [1, 4, 8];
    audit.CurrentGDPGrowth = currentGDPGrowth(audit.GDPGrowthHorizons);
    audit.LegacyGDPGrowth = legacyGDPGrowth(audit.GDPGrowthHorizons);
end

function localPrintLegacyAudit(audit)
    fprintf("\nLegacy GLP MAP benchmark audit:\n");
    fprintf("  %s\n", audit.Note);
    fprintf("  Relative coefficient difference: %.4g\n", ...
        audit.RelativeCoefficientDifference);
    fprintf("  Relative sigma difference:       %.4g\n", ...
        audit.RelativeSigmaDifference);
    fprintf("  Relative forecast difference:    %.4g\n", ...
        audit.RelativeForecastDifference);
    fprintf("  Hyperparameters:\n");
    fprintf("    lambda1 current / legacy lambda: %.6g / %.6g\n", ...
        audit.CurrentHyperparameters.lambda1, ...
        audit.LegacyHyperparameters.lambda);
    fprintf("    lambda4 current / legacy miu:    %.6g / %.6g\n", ...
        audit.CurrentHyperparameters.lambda4SumOfCoefficients, ...
        audit.LegacyHyperparameters.miuSumOfCoefficients);
    fprintf("    lambda5 current / legacy theta:  %.6g / %.6g\n", ...
        audit.CurrentHyperparameters.lambda5DummyInitialObservation, ...
        audit.LegacyHyperparameters.thetaDummyInitialObservation);
    fprintf("    alpha equivalent current/legacy: %.6g / %.6g\n", ...
        audit.CurrentHyperparameters.alphaEquivalent, ...
        audit.LegacyHyperparameters.alpha);
    fprintf("  GDP growth at mode, h = 1/4/8:\n");
    fprintf("    current: %s\n", mat2str(audit.CurrentGDPGrowth.', 5));
    fprintf("    legacy:  %s\n", mat2str(audit.LegacyGDPGrowth.', 5));
end

function legacyRun = localQuietLegacyGLP(Y, numLags)
    legacyFcn = @() bvarGLP(Y, numLags, 'mcmc', 0); %#ok<NASGU>
    [~, legacyRun] = evalc("legacyFcn()");
end

function r = localRelativeNorm(x, y)
    r = norm(x(:) - y(:)) / max(norm(y(:)), eps);
end

function localSaveFigure(fig, outputFile)
    exportgraphics(fig, outputFile, Resolution=150);
    fprintf("Saved figure to %s\n", outputFile);
end

function scaled = localScaleIRF(values, variable)
    if variable < 7
        scaled = values/4*100;
    else
        scaled = values*100;
    end
end

function label = localQuarterLabel(serialDate, yearDigits)
    dt = datetime(serialDate, ConvertFrom="datenum");
    if yearDigits == 2
        label = sprintf("Q%d-%02d", quarter(dt), mod(year(dt), 100));
    else
        label = sprintf("Q%d-%04d", quarter(dt), year(dt));
    end
end
