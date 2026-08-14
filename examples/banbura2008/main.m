%MAIN Replicate Banbura, Giannone and Reichlin (2008), "Large Bayesian VARs".
%
%   The script targets the paper's main empirical exercises:
%     * Tables 1 and 4: rolling relative MSFE for SMALL, CEE, MEDIUM, LARGE.
%     * Figure 1 and Table 5: recursive monetary-policy IRFs and FEVDs.
%
%   The paper PDF is read with Text Analytics Toolbox. The Stock-Watson data
%   are not shipped; place banbura2008Data.mat beside this file with a
%   monthly timetable named tt. See README.md for the expected variable names.
%
%   See also READBANBURA2008PAPER, BANBURA2008MODELSPECIFICATIONS,
%   MINNESOTABVARM, SVAR.IRF, SVAR.FEVD.

clear
clc

thisFolder = fileparts(mfilename("fullpath"));
repoRoot = fileparts(fileparts(thisFolder));

%% Paper metadata from Text Analytics Toolbox
paper = readBanbura2008Paper(fullfile(repoRoot, "RePEc_eca_wpaper_2008_033.pdf"));
fprintf("%s, %s\n", paper.Title, paper.WorkingPaper);
fprintf("Extracted %d characters from the PDF with Text Analytics Toolbox.\n\n", ...
    paper.NumCharacters);

%% Settings
fastMode = strcmpi(getenv("BGR_FAST"), "1");
numLags = paper.NumLags;
horizons = paper.ForecastHorizons;
numDraws = 1000;
targetVariables = ["EMPL" "CPI" "FFR"];
runLarge = true;
maxForecastOrigins = Inf;
resultsFile = fullfile(thisFolder, "banbura2008Results.mat");
rng(0, "twister");

if fastMode
    horizons = [1 12];
    numDraws = 100;
    runLarge = false;
    maxForecastOrigins = 12;
    fprintf("Fast mode enabled: reduced horizons, draws, and forecast origins.\n\n");
end

%% Data
dataFile = fullfile(thisFolder, "banbura2008Data.mat");
if ~isfile(dataFile)
    error("banbura2008:missingData", ...
        "Place banbura2008Data.mat in %s with a monthly timetable named tt.", ...
        thisFolder);
end

loaded = load(dataFile);
if ~isfield(loaded, "tt") || ~istimetable(loaded.tt)
    error("banbura2008:invalidData", ...
        "banbura2008Data.mat must contain a monthly timetable named tt.");
end
tt = sortrows(loaded.tt);

if isfield(loaded, "priorMean")
    priorMeanByName = localPriorMeanMap(tt, loaded.priorMean);
else
    priorMeanByName = dictionary(string(tt.Properties.VariableNames), ...
        ones(1, width(tt)));
end

%% Model panels
specs = banbura2008ModelSpecifications();
modelNames = ["SMALL" "CEE" "MEDIUM"];
if runLarge && width(tt) >= specs.LARGE.RequiredVariableCount
    modelNames(end+1) = "LARGE";
end

%% Forecast evaluation: Tables 1 and 4
fprintf("Rolling forecast evaluation, relative to random walk with drift.\n");
fprintf("Evaluation sample: %s to %s. Lags: %d. Window: 120 months.\n\n", ...
    string(paper.ForecastEvaluationStart, "yyyy-MM"), ...
    string(paper.ForecastEvaluationEnd, "yyyy-MM"), numLags);

withoutSOC = localForecastTable(tt, specs, modelNames, targetVariables, ...
    horizons, numLags, paper.ForecastEvaluationStart, ...
    paper.ForecastEvaluationEnd, priorMeanByName, ...
    IncludeSumOfCoefficients=false, MaxForecastOrigins=maxForecastOrigins);
withSOC = localForecastTable(tt, specs, modelNames, targetVariables, ...
    horizons, numLags, paper.ForecastEvaluationStart, ...
    paper.ForecastEvaluationEnd, priorMeanByName, ...
    IncludeSumOfCoefficients=true, MaxForecastOrigins=maxForecastOrigins);

fprintf("Table 1 style: Minnesota BVAR only\n");
disp(withoutSOC);
fprintf("Table 4 style: with sum-of-coefficients prior, lambda4 = 10*lambda1\n");
disp(withSOC);

%% Structural analysis: Figure 1 and Table 5
structuralTT = tt(timerange(paper.StructuralSampleStart, ...
    paper.StructuralSampleEnd, "closed"), :);

for modelName = modelNames
    spec = localResolveSpec(specs.(modelName), structuralTT);
    Y = structuralTT{:, spec.Variables};
    priorMean = localPriorMean(priorMeanByName, spec.Variables);

    posterior = localEstimateBVAR(Y, spec, numLags, priorMean, ...
        IncludeSumOfCoefficients=true);
    pointVAR = bvar2var(posterior);

    % Cholesky ordering: slow-moving variables, FFR, fast-moving variables.
    impact = chol(pointVAR.Covariance, "lower");
    scaledImpact = impact(:, spec.PolicyIndex) ./ impact(spec.PolicyIndex, spec.PolicyIndex);

    pointIRF = svar.irf(pointVAR, scaledImpact, max(paper.StructuralHorizons));
    pointFEVD = svar.fevd(pointVAR, impact, max(paper.StructuralHorizons));

    [coeffDraws, sigmaDraws] = simulate(posterior, NumDraws=numDraws);
    irfDraws = zeros(numel(paper.StructuralHorizons), numel(spec.Variables), numDraws);
    for draw = 1:numDraws
        drawVAR = svar.varmFromCoefficients(posterior, ...
            reshape(coeffDraws(:,draw), [], numel(spec.Variables)), sigmaDraws(:,:,draw));
        drawImpact = chol(sigmaDraws(:,:,draw), "lower");
        drawShock = drawImpact(:, spec.PolicyIndex) ...
            ./ drawImpact(spec.PolicyIndex, spec.PolicyIndex);
        irfDraws(:,:,draw) = svar.irf(drawVAR, drawShock, max(paper.StructuralHorizons));
    end

    localPlotPolicyIRF(modelName, paper.StructuralHorizons, pointIRF, irfDraws, ...
        spec.Variables, targetVariables);
    localPrintFEVD(modelName, pointFEVD, spec, targetVariables, ...
        [1 3 6 12 24 36 48]);
end

save(resultsFile, "paper", "withoutSOC", "withSOC", "modelNames", ...
    "horizons", "numLags", "fastMode");
fprintf("\nSaved replication tables to %s\n", resultsFile);

%% Local functions
function results = localForecastTable(tt, specs, modelNames, targetVariables, ...
        horizons, numLags, evalStart, evalEnd, priorMeanByName, nvp)
arguments
    tt timetable
    specs struct
    modelNames (1,:) string
    targetVariables (1,:) string
    horizons (1,:) double
    numLags (1,1) double
    evalStart (1,1) datetime
    evalEnd (1,1) datetime
    priorMeanByName dictionary
    nvp.IncludeSumOfCoefficients (1,1) logical = false
    nvp.MaxForecastOrigins (1,1) double {mustBePositive} = Inf
end

rows = strings(0,1);
modelColumn = strings(0,1);
hColumn = zeros(0,1);
relativeMSFE = zeros(0,1);

for modelName = modelNames
    spec = localResolveSpec(specs.(modelName), tt);
    modelTT = tt(:, spec.Variables);
    Y = modelTT{:,:};
    dates = modelTT.Time;
    targetIndex = arrayfun(@(v) find(spec.Variables == v, 1), targetVariables);
    priorMean = localPriorMean(priorMeanByName, spec.Variables);

    for h = horizons
        forecastErrors = zeros(0, numel(targetVariables));
        benchmarkErrors = zeros(0, numel(targetVariables));
        lastOrigin = evalEnd - calmonths(h);
        origins = dates(dates >= evalStart - calmonths(12-h) & dates <= lastOrigin);
        if isfinite(nvp.MaxForecastOrigins) && numel(origins) > nvp.MaxForecastOrigins
            origins = origins(1:nvp.MaxForecastOrigins);
        end

        for origin = reshape(origins, 1, [])
            originIndex = find(dates == origin, 1);
            if originIndex < 120 || originIndex + h > height(modelTT)
                continue
            end

            window = Y(originIndex-119:originIndex, :);
            actual = Y(originIndex+h, targetIndex);

            posterior = localEstimateBVAR(window, spec, numLags, priorMean, ...
                IncludeSumOfCoefficients=nvp.IncludeSumOfCoefficients);
            pointVAR = bvar2var(posterior);
            forecastPath = forecast(pointVAR, h, window);
            predicted = forecastPath(end, targetIndex);
            benchmark = localRandomWalkWithDrift(window(:, targetIndex), h);

            forecastErrors(end+1,:) = predicted - actual; %#ok<AGROW>
            benchmarkErrors(end+1,:) = benchmark - actual; %#ok<AGROW>
        end

        msfe = mean(forecastErrors.^2, 1, "omitnan");
        benchmarkMSFE = mean(benchmarkErrors.^2, 1, "omitnan");
        ratio = msfe ./ benchmarkMSFE;

        rows = [rows; targetVariables(:)]; %#ok<AGROW>
        modelColumn = [modelColumn; repmat(modelName, numel(targetVariables), 1)]; %#ok<AGROW>
        hColumn = [hColumn; repmat(h, numel(targetVariables), 1)]; %#ok<AGROW>
        relativeMSFE = [relativeMSFE; ratio(:)]; %#ok<AGROW>
    end
end

results = table(modelColumn, hColumn, rows, relativeMSFE, ...
    VariableNames=["Model" "Horizon" "Variable" "RelativeMSFE"]);
end

function posterior = localEstimateBVAR(Y, spec, numLags, priorMean, nvp)
arguments
    Y double
    spec struct
    numLags (1,1) double
    priorMean (1,:) double
    nvp.IncludeSumOfCoefficients (1,1) logical = false
end

numSeries = size(Y, 2);
if isinf(spec.lambda1)
    prior = diffusebvarm(numSeries, numLags, SeriesNames=spec.Variables);
    posterior = estimate(prior, Y, Display="off");
    return
end

args = {"lambda1", spec.lambda1, "lambda3", 1, "PriorMean", priorMean, ...
    "SeriesNames", spec.Variables, "Psi", "conditional"};
if nvp.IncludeSumOfCoefficients
    args = [args, {"lambda4", spec.lambda4}];
end
prior = minnesotabvarm(numSeries, numLags, Y, args{:});
posterior = estimate(prior, Display="off");
end

function spec = localResolveSpec(spec, tt)
if spec.Name == "LARGE"
    names = string(tt.Properties.VariableNames);
    if ~any(names == spec.PolicyVariable)
        error("banbura2008:missingPolicyVariable", ...
            "The LARGE data must include %s.", spec.PolicyVariable);
    end
    slow = names([1:70 108:131]);
    fast = names([71:85 87:107]);
    spec.Variables = [slow spec.PolicyVariable fast];
    spec.PolicyIndex = numel(slow) + 1;
end

missing = setdiff(spec.Variables, string(tt.Properties.VariableNames));
if ~isempty(missing)
    error("banbura2008:missingVariables", ...
        "Missing variables for %s: %s", spec.Name, strjoin(missing, ", "));
end
end

function map = localPriorMeanMap(tt, priorMean)
names = string(tt.Properties.VariableNames);
if numel(priorMean) ~= numel(names)
    error("banbura2008:invalidPriorMean", ...
        "priorMean must have one value per timetable variable.");
end
map = dictionary(names, reshape(double(priorMean), 1, []));
end

function priorMean = localPriorMean(map, names)
priorMean = ones(1, numel(names));
for j = 1:numel(names)
    if isKey(map, names(j))
        priorMean(j) = map(names(j));
    end
end
end

function forecast = localRandomWalkWithDrift(Y, h)
drift = mean(diff(Y, 1, 1), 1, "omitnan");
forecast = Y(end,:) + h*drift;
end

function localPlotPolicyIRF(modelName, horizons, pointIRF, irfDraws, names, targetVariables)
targetIndex = arrayfun(@(v) find(names == v, 1), targetVariables);
band = quantile(irfDraws(:, targetIndex, :), [0.05 0.16 0.84 0.95], 3);

figure(Color="w");
layout = tiledlayout(numel(targetVariables), 1, TileSpacing="compact", Padding="compact");
for j = 1:numel(targetVariables)
    ax = nexttile(layout);
    fill(ax, [horizons fliplr(horizons)], ...
        [band(:,j,1); flipud(band(:,j,4))]', [0.88 0.90 0.95], ...
        EdgeColor="none");
    hold(ax, "on");
    fill(ax, [horizons fliplr(horizons)], ...
        [band(:,j,2); flipud(band(:,j,3))]', [0.74 0.80 0.92], ...
        EdgeColor="none");
    plot(ax, horizons, pointIRF(:, targetIndex(j)), "k", LineWidth=1.2);
    yline(ax, 0, ":");
    xlim(ax, [horizons(1) horizons(end)]);
    title(ax, targetVariables(j));
    box(ax, "on");
end
title(layout, modelName + ": response to a 100 bp federal funds rate shock");
end

function localPrintFEVD(modelName, fevd, spec, targetVariables, horizons)
fprintf("\n%s FEVD share of the monetary-policy shock, percent\n", modelName);
for target = targetVariables
    targetIndex = find(spec.Variables == target, 1);
    shares = 100*fevd(horizons + 1, targetIndex, spec.PolicyIndex);
    fprintf("  %-5s %s\n", target, mat2str(round(shares(:)', 1)));
end
end
