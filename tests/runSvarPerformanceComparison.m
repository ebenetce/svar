function [summary, results, plotFiles] = runSvarPerformanceComparison(opts)
%RUNSVARPERFORMANCECOMPARISON Compare SVAR performance against built-ins.
%   SUMMARY = RUNSVARPERFORMANCECOMPARISON runs the SVAR IRF and FEVD
%   performance tests, checks that the SVAR implementations are faster than
%   orthogonalized and generalized built-in baselines, and writes comparison
%   plots to public/performance.
%
%   [SUMMARY,RESULTS,PLOTFILES] = RUNSVARPERFORMANCECOMPARISON returns the
%   comparison table, raw TimeResult array, and generated plot file names.

arguments
    opts.SampleSize (1,1) double {mustBeInteger, mustBePositive} = 6
    opts.MinimumSpeedup (1,1) double {mustBePositive} = 1
    opts.DoPlots (1,1) logical = false
end

root = currentProject().RootFolder;
OutputFolder = fullfile(root, "public", "performance");

if ~isfolder(OutputFolder)
    mkdir(OutputFolder);
end

import matlab.perftest.TimeExperiment
suite = [testsuite(fullfile(root, "tests", "svarIrfPerformanceTest.m")), ...
    testsuite(fullfile(root, "tests", "svarFevdPerformanceTest.m"))];
experiment = TimeExperiment.withFixedSampleSize(opts.SampleSize);
results = run(experiment, suite);

summary = localSummary(results);
if any(summary.Speedup < opts.MinimumSpeedup)
    disp(summary);
    error("SVAR:PerformanceComparison:NotFaster", ...
        "SVAR median timing did not beat every built-in baseline.");
end

if opts.DoPlots
    plotFiles = localComparisonPlots(results, OutputFolder);
end
disp(summary);

end

function summary = localSummary(results)
pairs = [
    "IRF"  "armairf orthogonalized"  "testSvarIrf"   "testArmaIrfOrthogonalized"
    "IRF"  "armairf generalized"     "testSvarIrf"   "testArmaIrfGeneralized"
    "IRF"  "varm/irf orthogonalized" "testSvarIrf"   "testVarmIrfOrthogonalized"
    "IRF"  "varm/irf generalized"    "testSvarIrf"   "testVarmIrfGeneralized"
    "FEVD" "armafevd orthogonalized" "testSvarFevd"  "testArmaFevdOrthogonalized"
    "FEVD" "armafevd generalized"    "testSvarFevd"  "testArmaFevdGeneralized"
    "FEVD" "varm/fevd orthogonalized" "testSvarFevd" "testVarmFevdOrthogonalized"
    "FEVD" "varm/fevd generalized"   "testSvarFevd"  "testVarmFevdGeneralized"
    ];

numPairs = size(pairs, 1);
Family = strings(numPairs, 1);
Baseline = strings(numPairs, 1);
SvarMedianSeconds = zeros(numPairs, 1);
BaselineMedianSeconds = zeros(numPairs, 1);
Speedup = zeros(numPairs, 1);

for pairIndex = 1:numPairs
    svarResult = localResultByMethod(results, pairs(pairIndex, 3));
    baselineResult = localResultByMethod(results, pairs(pairIndex, 4));
    svarMedian = median(svarResult.Samples.MeasuredTime);
    baselineMedian = median(baselineResult.Samples.MeasuredTime);

    Family(pairIndex) = pairs(pairIndex, 1);
    Baseline(pairIndex) = pairs(pairIndex, 2);
    SvarMedianSeconds(pairIndex) = svarMedian;
    BaselineMedianSeconds(pairIndex) = baselineMedian;
    Speedup(pairIndex) = baselineMedian/svarMedian;
end

summary = table(Family, Baseline, SvarMedianSeconds, ...
    BaselineMedianSeconds, Speedup);
end

function plotFiles = localComparisonPlots(results, outputFolder)
families = ["IRF" "FEVD"];
plotFiles = strings(size(families));

for familyIndex = 1:numel(families)
    family = families(familyIndex);
    fig = figure();

    switch family
        case "IRF"
            svarMethod = "testSvarIrf";
            baselineMethods = ["testArmaIrfOrthogonalized", ...
                "testArmaIrfGeneralized", "testVarmIrfOrthogonalized", ...
                "testVarmIrfGeneralized"];
        case "FEVD"
            svarMethod = "testSvarFevd";
            baselineMethods = ["testArmaFevdOrthogonalized", ...
                "testArmaFevdGeneralized", "testVarmFevdOrthogonalized", ...
                "testVarmFevdGeneralized"];
    end

    svarResult = localResultByMethod(results, svarMethod);
    baselineResults = arrayfun(@(methodName) ...
        localResultByMethod(results, methodName), baselineMethods);
    measurementResults = repmat(svarResult, size(baselineResults));
    comparisonPlot(baselineResults, measurementResults, "median", ...
        Parent=fig, Scale="linear");
    title(family + " performance: built-in baseline vs SVAR");

    plotFiles(familyIndex) = fullfile(outputFolder, ...
        lower(family) + "-performance-comparison.png");
end
end

function result = localResultByMethod(results, methodName)
names = string({results.Name});
matches = contains(names, methodName);
if nnz(matches) ~= 1
    error("SVAR:PerformanceComparison:ResultNotFound", ...
        "Expected one performance result matching '%s', found %d.", ...
        methodName, nnz(matches));
end
result = results(matches);
end
