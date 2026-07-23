%MAIN Replicate Giannone, Lenza and Primiceri (2012) with the SVAR toolbox.
%
%   "Prior Selection for Vector Autoregressions", NBER WP 18467 (GLP_2012.pdf).
%
%   The paper treats the tightness of a conjugate Minnesota / sum-of-
%   coefficients / dummy-initial-observation prior as unknown and infers it
%   from the marginal likelihood. This script reproduces that exercise on the
%   authors' own MEDIUM dataset (Stock-Watson 2008, 7 series, 1959Q1-2008Q4)
%   using GLP, SVAR.MINNESOTAMNIWBVARM and the Econometrics Toolbox, and
%   audits the result against the authors' original code in lenzaPrimiceri/.
%
%   What is targeted
%   ----------------
%     * Figure 1 of the paper - the posterior of lambda in a Minnesota-only
%       BVAR, for the SMALL and MEDIUM datasets.
%     * The authors' released scripts ExamplePredictiveDensity.m and
%       ExampleIRFs.m, which run the full three-prior model on the MEDIUM
%       dataset. These are the concrete numeric target; the paper's own
%       figures 2-5 use the 22-variable LARGE dataset, which is not shipped.
%
%   Hyperparameter mapping (paper -> toolbox)
%   -----------------------------------------
%       lambda  overall Minnesota tightness         -> lambda1
%       alpha   lag decay, variance ~ 1/l^alpha     -> lambda3 = alpha/2
%       mu      sum-of-coefficients tightness       -> lambda4
%       delta   dummy-initial-observation tightness -> lambda5
%       psi     Inverse-Wishart scale diagonal      -> Psi
%       d       Inverse-Wishart d.o.f., fixed n+2   -> DoF (set by the class)
%
%   Hyperparameter uncertainty is integrated out, not conditioned away: the
%   GLP call below runs the Metropolis step of the paper's appendix B, so the
%   bands reported here reflect uncertainty about lambda, mu, delta and psi as
%   well as about the VAR parameters. The script also reports what the bands
%   would have been WITHOUT that, which is the paper's whole point.
%
%   Runtime is roughly 10 minutes, dominated by the two per-draw loops.
%   See the "Replication gaps" section at the end for what still differs.
%
%   See also GLP, SVAR.MINNESOTAMNIWBVARM, HYPERPRIOR, LOGMARGINALLIKELIHOOD.

clear
clc

thisFolder = fileparts(mfilename("fullpath"));

%% Settings
numLags          = 5;                % paper: 5 lags in every VAR
forecastHorizons = 1:8;              % quarters ahead
irfHorizon       = 20;               % IRF horizons plotted (0:19)
policyShock      = 7;                % FedFunds, ordered last -> Cholesky MP shock
numDraws         = 10000;            % GLP keep 10000 of 20000 Metropolis draws
burnIn           = 10000;            % ... discarding the first half
proposalScale    = 0.8;              % tuned for ~25% acceptance (see below)
randomSeed       = 0;
quantiles        = [0.16 0.50 0.84]; % bands reported by the authors' code
runLegacyMCMC    = false;            % compare chains with bvarGLP (adds ~2 min)

%% Data
% Annualised log-levels (4 x log) for the six real/nominal series, the funds
% rate in levels divided by 100.
load(fullfile(thisFolder, "data.mat"), "tt");

Y           = tt{:,:};
seriesNames = string(tt.Properties.VariableNames);
numSeries   = width(tt);
lastDate    = tt.Time(end);

%% Hyperpriors
% Gamma densities with the modes recommended by Sims and Zha (1998) and the
% standard deviations in section 3 of the paper.
lambda1Prior = hyperprior("Gamma", 0.2, 0.4, Bounds = [1e-4, 5]);   % lambda
lambda4Prior = hyperprior("Gamma", 1,   1,   Bounds = [1e-4, 50]);  % mu
lambda5Prior = hyperprior("Gamma", 1,   1,   Bounds = [1e-4, 50]);  % delta

% psi carries an Inverse-Gamma with shape = scale = 0.02^2, calibrated in the
% paper for data in 4 x logs. That density has neither a mean nor a variance,
% so it is unreachable from the (mode, sd) parameterisation and has to be
% given natively. Each series is started at its AR(1) residual variance and
% searched over [psi/100, 100*psi], the authors' bounds.
psi0     = svar.estimateResidualVariances(Y, 1, Method = "conditional");
psiPrior = arrayfun(@(v) hyperprior("InverseGamma", 0.02^2, 0.02^2, ...
    Parameterization = "native", X0 = v, Bounds = [v/100, v*100]), psi0);

optimOptions = optimoptions("fmincon", ...
    Display                  = "off", ...
    FiniteDifferenceStepSize = 1e-4, ...
    FunctionTolerance        = 1e-12, ...
    StepTolerance            = 1e-12, ...
    ConstraintTolerance      = 1e-12, ...
    MaxIterations            = 1000, ...
    MaxFunctionEvaluations   = 20000);

%% Figure 1: posterior of lambda in a Minnesota-only BVAR
% Section 4.3 simplifies the model to a Minnesota prior alone, with psi fixed
% at the AR(1) residual variances. That leaves a SINGLE free hyperparameter,
% so its posterior can be normalised on a grid - no MCMC needed - and the
% mode/spread can be compared with the paper's figure 1 directly.
smallModel = [1 2 7];   % GDP, GDP deflator, federal funds rate

lambdaGrid = linspace(0.02, 1.6, 200);
smallPost  = localLambdaPosterior(Y(:, smallModel), numLags, lambdaGrid, lambda1Prior);
mediumPost = localLambdaPosterior(Y, numLags, lambdaGrid, lambda1Prior);

figure(Color = "w");
plot(lambdaGrid, smallPost, LineWidth = 1.4, DisplayName = "SMALL (3 variables)");
hold on
plot(lambdaGrid, mediumPost, LineWidth = 1.4, DisplayName = "MEDIUM (7 variables)");
plot(lambdaGrid, exp(lambda1Prior.logpdf(lambdaGrid)), "k--", ...
    LineWidth = 1.1, DisplayName = "Hyperprior");
xlabel("\lambda"); ylabel("Density");
title("Posterior of the Minnesota tightness (paper figure 1)");
legend(Location = "northeast"); box on

fprintf("Posterior of lambda, Minnesota-only model (paper figure 1):\n");
fprintf("                     mode    peak density    paper figure 1\n");
fprintf("  SMALL  (3 series) %6.3f %13.1f    mode ~0.42, peak ~5.3\n", ...
    localGridMode(lambdaGrid, smallPost), max(smallPost));
fprintf("  MEDIUM (7 series) %6.3f %13.1f    mode ~0.17, peak ~19\n", ...
    localGridMode(lambdaGrid, mediumPost), max(mediumPost));
fprintf("The mode and the spread both shrink with the size of the model.\n\n");

%% Hyperparameter posterior, full model
% All three priors active: Minnesota, sum-of-coefficients, dummy initial
% observation. lambda3 is FIXED at 1, i.e. alpha = 2 (the paper's default).
%
% NumDraws switches on the paper's appendix-B Metropolis step: the maximiser
% locates the mode and supplies the proposal covariance, then the chain
% explores around it. Ten free hyperparameters here (lambda, mu, delta and
% seven psi), so the classic 2.38/sqrt(d) ~ 0.75 starting point for the step
% size is about right.
rng(randomSeed, "twister");
[priorMdl, searchInfo, chain] = glp(numSeries, numLags, tt, ...
    lambda1       = lambda1Prior, ...
    lambda3       = 1, ...
    lambda4       = lambda4Prior, ...
    lambda5       = lambda5Prior, ...
    Vc            = 10e6, ...
    Psi           = psiPrior, ...
    OptimOptions  = optimOptions, ...
    NumDraws      = numDraws, ...
    BurnIn        = burnIn, ...
    ProposalScale = proposalScale);

fprintf("Hyperparameter posterior mode (all priors):\n");
fprintf("  lambda (lambda1) = %.6f\n", priorMdl.lambda1);
fprintf("  mu     (lambda4) = %.6f\n", priorMdl.lambda4);
fprintf("  delta  (lambda5) = %.6f\n", priorMdl.lambda5);
fprintf("  alpha  (2*l3)    = %.6f\n", 2*priorMdl.lambda3);
fprintf("  psi              = %s\n", mat2str(priorMdl.ResidualVariances, 5));
fprintf("  log posterior    = %.4f (exitflag %d)\n\n", ...
    -searchInfo.Objective, searchInfo.ExitFlag);

% The mode is not the whole story: these are the marginals the paper's
% hierarchical argument is about.
hyperTable = table( ...
    [priorMdl.lambda1; priorMdl.lambda4; priorMdl.lambda5], ...
    [mean(chain.lambda1); mean(chain.lambda4); mean(chain.lambda5)], ...
    [std(chain.lambda1);  std(chain.lambda4);  std(chain.lambda5)], ...
    [quantile(chain.lambda1, 0.16); quantile(chain.lambda4, 0.16); quantile(chain.lambda5, 0.16)], ...
    [quantile(chain.lambda1, 0.84); quantile(chain.lambda4, 0.84); quantile(chain.lambda5, 0.84)], ...
    VariableNames = ["Mode", "Mean", "Std", "Q16", "Q84"], ...
    RowNames = ["lambda", "mu (soc)", "delta (dio)"]);

fprintf("Hyperparameter posterior (%d draws, acceptance %.3f):\n", ...
    chain.NumDraws, chain.AcceptanceRate);
disp(hyperTable);

%% Posterior of the VAR parameters at the modal hyperparameters
posterior = estimate(priorMdl, Display = "off");

modeCoefficients = reshape(posterior.Mu, [], numSeries);
modeSigma        = posterior.Omega/(posterior.DoF + numSeries + 1);  % IW mode
modeVAR          = svar.varmFromCoefficients(posterior, modeCoefficients, modeSigma);

%% Draws of the VAR parameters and of the predictive density
% chain.Coefficients / chain.Sigma already integrate over the hyperparameters:
% each draw carries its own lambda, mu, delta and psi.
presample = Y(end-numLags+1:end, :);

[irfDraws, forecastDraws] = localDrawPaths(posterior, chain.Coefficients, ...
    chain.Sigma, policyShock, irfHorizon, forecastHorizons, presample);

% The same exercise holding the hyperparameters AT the mode - what the model
% would have said before appendix B. Kept only to measure what integrating
% them out is worth.
[fixedCoefficients, fixedCovariances] = simulate(posterior, NumDraws = numDraws);
fixedIRFDraws = localDrawPaths(posterior, ...
    reshape(fixedCoefficients, [], numSeries, numDraws), fixedCovariances, ...
    policyShock, irfHorizon, forecastHorizons, presample);

%% Forecasts of average GDP growth (ExamplePredictiveDensity)
% Annualised average growth over h quarters: (y_{T+h} - y_T)/h, in percent.
modeForecast = forecast(modeVAR, max(forecastHorizons), presample);
modeGrowth   = localAverageGrowth(modeForecast(forecastHorizons, 1), Y(end,1), forecastHorizons);
drawGrowth   = localAverageGrowth(squeeze(forecastDraws(:,1,:)), Y(end,1), forecastHorizons);

reported = [1 4 8];
growthTable = table(reported.', modeGrowth(reported), ...
    median(drawGrowth(reported,:), 2), std(drawGrowth(reported,:), 0, 2), ...
    VariableNames = ["QuartersAhead", "PointAtMode", "MedianPredictive", "StdPredictive"]);

fprintf("Forecasting in %s (GDP annualised average growth)\n", ...
    string(lastDate, "QQQ-yy"));
disp(growthTable);

figure(Color = "w");
layout = tiledlayout(2, 1, TileSpacing = "compact", Padding = "compact");
for horizon = [1 4]
    ax = nexttile(layout);
    histogram(ax, drawGrowth(horizon,:), 50, EdgeColor = "none");
    xline(ax, modeGrowth(horizon), LineWidth = 1.25);
    title(ax, sprintf("Avg GDP growth %d qrts ahead (mode: %.2f)", ...
        horizon, modeGrowth(horizon)));
    xlabel(ax, "Annualised average growth, percent");
    xlim(ax, [-15 10]); box(ax, "on");
end
title(layout, "Predictive density, forecasting in " + string(lastDate, "QQQ-yy"));

%% Impulse responses to a monetary policy shock (ExampleIRFs)
% Cholesky identification with the funds rate ordered last.
modeIRF = svar.irf(modeVAR, chol(modeSigma, "lower"), irfHorizon - 1);
modeIRF = localScaleIRF(modeIRF(:,:,policyShock));
bandIRF = localScaleIRF(quantile(irfDraws, quantiles, 3));

localPlotIRFs(modeIRF, bandIRF, seriesNames, quantiles, ...
    "Monetary policy shock, all priors");

fprintf("\nImpulse responses to a 1 s.d. monetary policy tightening:\n");
fprintf("  Funds rate on impact       : %+7.2f percentage points\n", modeIRF(1, policyShock));
fprintf("  GDP, avg over horizons 1-8 : %+7.2f percent\n", mean(modeIRF(2:9, 1)));
fprintf("  Deflator, avg horizons 4-12: %+7.2f percent\n", mean(modeIRF(5:13, 2)));

%% What integrating over the hyperparameters is worth
% Section 2 of the paper: conditioning on a single value of gamma understates
% uncertainty. This is that claim, measured on this dataset.
fixedBand = localScaleIRF(quantile(fixedIRFDraws, quantiles, 3));
widthRatio = (bandIRF(:,:,3) - bandIRF(:,:,1)) ...
    ./ (fixedBand(:,:,3) - fixedBand(:,:,1));

% Horizon 0 is degenerate for every variable ordered before the shock: a
% lower-triangular impact matrix pins those responses at exactly zero, so both
% band widths are zero and the ratio is 0/0.
widthRatio(1, 1:policyShock-1) = NaN;

fprintf("\n68%% IRF band width, integrated over hyperparameters vs at the mode\n");
fprintf("(ratio > 1 means conditioning on the modal hyperparameters is overconfident)\n");
disp(array2table(median(widthRatio, 1, "omitnan"), VariableNames = seriesNames, ...
    RowNames = "median ratio over horizons"));

%% Audit against the authors' original code
% bvarGLP is the authors' implementation. Both sides now carry the same
% hyperpriors, so the log posteriors are directly comparable; the two remaining
% differences are the definition of ybar (see the gaps section) and csminwel
% versus fmincon.
oldPath = path;
cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(thisFolder, "lenzaPrimiceri"));

[~, legacy] = evalc("bvarGLP(Y, numLags, 'mcmc', 0)");
legacyMode  = legacy.postmax;

% The authors order the constant FIRST; the toolbox orders it after the lags.
toolboxInLegacyOrder = [modeCoefficients(end,:); modeCoefficients(1:end-1,:)];
legacyGrowth = localAverageGrowth(legacyMode.forecast(:,1), Y(end,1), forecastHorizons);

fprintf("\nAudit against bvarGLP (authors' code)\n");
fprintf("                        toolbox      authors\n");
fprintf("  lambda            %10.6f   %10.6f\n", priorMdl.lambda1, legacyMode.lambda);
fprintf("  mu                %10.6f   %10.6f\n", priorMdl.lambda4, legacyMode.miu);
fprintf("  delta             %10.6f   %10.6f\n", priorMdl.lambda5, legacyMode.theta);
fprintf("  alpha             %10.6f   %10.6f\n", 2*priorMdl.lambda3, legacyMode.alpha);
fprintf("  log posterior     %10.4f   %10.4f\n", -searchInfo.Objective, legacyMode.logPost);
fprintf("  relative difference in coefficients : %.3g\n", ...
    localRelativeNorm(toolboxInLegacyOrder, legacyMode.betahat));
fprintf("  relative difference in Sigma        : %.3g\n", ...
    localRelativeNorm(modeSigma, legacyMode.sigmahat));
fprintf("  relative difference in forecasts    : %.3g\n", ...
    localRelativeNorm(modeForecast(forecastHorizons,:), legacyMode.forecast));
fprintf("  GDP growth at mode, h = 1/4/8, toolbox: %s\n", mat2str(modeGrowth(reported).', 5));
fprintf("  GDP growth at mode, h = 1/4/8, authors: %s\n", mat2str(legacyGrowth(reported).', 5));

if runLegacyMCMC
    % Same comparison for the Metropolis step itself. MCMCconst is the
    % authors' recommended value for this model.
    [~, legacyChain] = evalc("bvarGLP(Y, numLags, 'mcmc', 1, 'MCMCconst', 1.6)"); %#ok<UNRCH>
    legacyChain = legacyChain.mcmc;

    fprintf("\n  Hyperparameter chains        toolbox      authors\n");
    fprintf("    acceptance rate       %10.3f   %10.3f\n", ...
        chain.AcceptanceRate, legacyChain.ACCrate);
    fprintf("    distinct lambda draws %10d   %10d\n", ...
        numel(unique(chain.lambda1)), numel(unique(legacyChain.lambda)));
    fprintf("    lambda  mean / sd     %5.3f / %-5.3f %5.3f / %.3f\n", ...
        mean(chain.lambda1), std(chain.lambda1), ...
        mean(legacyChain.lambda), std(legacyChain.lambda));
    fprintf("    mu      mean / sd     %5.3f / %-5.3f %5.3f / %.3f\n", ...
        mean(chain.lambda4), std(chain.lambda4), ...
        mean(legacyChain.miu), std(legacyChain.miu));
    fprintf("    delta   mean / sd     %5.3f / %-5.3f %5.3f / %.3f\n", ...
        mean(chain.lambda5), std(chain.lambda5), ...
        mean(legacyChain.theta), std(legacyChain.theta));
end

%% Replication gaps
% What GLP (2012) does that the toolbox still cannot reproduce.
%
% 1. THE pos OPTION. bvarGLP takes `pos`, the variables entering in first
%    differences, and zeroes BOTH their Minnesota prior mean on the own first
%    lag AND the corresponding block of the sum-of-coefficients dummy
%    (ydnoc(pos,pos) = 0). PriorMean covers the first half; DUMMYMATRICES has
%    no equivalent for the second. It does not bite here (pos is empty in
%    every GLP example) but blocks any replication using differenced data.
%
% 2. ybar IN THE DUMMY OBSERVATIONS - a difference, not a deficiency. The
%    paper (p.9) defines ybar as "the average of the first p observations",
%    which is what SVAR.MINNESOTAMNIWBVARM does; the authors' code trims the
%    presample first and so averages observations p+1 ... 2p. The toolbox is
%    faithful to the paper and the released code is not. It accounts for
%    0.085 of the log-posterior difference reported by the audit above.
%
% Everything else matches. Given the same hyperparameters, the toolbox log
% marginal likelihood agrees with the authors' formula to machine precision
% once ybar is aligned - the substantive verification of this port.
%
% On the Metropolis step, one difference is worth recording because it favours
% this implementation. The authors maximise in logit-transformed coordinates,
% so their Hessian has to be pushed back through a Jacobian before it can be
% used as a proposal covariance; GLP maximises in natural coordinates and
% needs no such step. Their recommended MCMCconst = 1.6 delivers a 3 per cent
% acceptance rate on this dataset - far below the 20-25 per cent their own
% paper and setpriors.m target - and only about 300 distinct hyperparameter
% values in 10000 draws. Set runLegacyMCMC above to see both chains side by
% side; the posterior moments agree, but the chains do not mix equally well.

%% Local functions
function posteriorDensity = localLambdaPosterior(Y, numLags, lambdaGrid, lambdaPrior)
%LOCALLAMBDAPOSTERIOR Normalised posterior of lambda in a Minnesota-only BVAR.
%   psi is FIXED at the AR(1) residual variances and the sum-of-coefficients
%   and dummy-initial-observation priors are off (lambda4 = lambda5 = Inf), so
%   lambda is the only unknown hyperparameter and its posterior is a
%   one-dimensional normalisation rather than an MCMC problem.
numSeries = size(Y, 2);
psi       = svar.estimateResidualVariances(Y, 1, Method = "conditional");

logPosterior = arrayfun(@(lambda) ...
    logMarginalLikelihood(svar.minnesotamniwbvarm(numSeries, numLags, Y, ...
        Psi = psi, lambda1 = lambda, lambda3 = 1, Vc = 10e6)) ...
    + lambdaPrior.logpdf(lambda), lambdaGrid);

density          = exp(logPosterior - max(logPosterior));
posteriorDensity = density/trapz(lambdaGrid, density);
end

function [irfDraws, forecastDraws] = localDrawPaths(template, coefficients, ...
    covariances, shock, irfHorizon, horizons, presample)
%LOCALDRAWPATHS Per-draw impulse responses and predictive forecast paths.
%   COEFFICIENTS is m-by-n-by-numDraws, the layout CHAIN returns. Each draw
%   becomes a VARM, which the Econometrics Toolbox then simulates forward;
%   the shock is identified by a Cholesky factor with the policy rate last.
numDraws  = size(covariances, 3);
numSeries = size(covariances, 1);

irfDraws      = zeros(irfHorizon, numSeries, numDraws);
forecastDraws = zeros(numel(horizons), numSeries, numDraws);

for draw = 1:numDraws
    sigma   = covariances(:,:,draw);
    impact  = chol(sigma, "lower");
    drawVAR = svar.varmFromCoefficients(template, coefficients(:,:,draw), sigma);

    irfDraws(:,:,draw)      = svar.irf(drawVAR, impact(:, shock), irfHorizon - 1);
    forecastDraws(:,:,draw) = simulate(drawVAR, max(horizons), Y0 = presample);
end
end

function modeValue = localGridMode(grid, density)
[~, at]   = max(density);
modeValue = grid(at);
end

function growth = localAverageGrowth(levelForecast, lastLevel, horizons)
%LOCALAVERAGEGROWTH Annualised average growth over h quarters, in percent.
growth = (levelForecast - lastLevel).*(100./horizons(:));
end

function scaled = localScaleIRF(irf)
%LOCALSCALEIRF Put every response in percent.
%   Series 1-6 are 4 x log, so a unit response is a quarter of a log point;
%   the funds rate is already a rate divided by 100.
scaled            = irf*100;
scaled(:,1:6,:)   = scaled(:,1:6,:)/4;
end

function localPlotIRFs(modeIRF, bandIRF, seriesNames, quantiles, plotTitle)
horizons = 0:(size(modeIRF, 1) - 1);
figure(Color = "w");
layout = tiledlayout(2, 4, TileSpacing = "compact", Padding = "compact");

for series = 1:numel(seriesNames)
    ax = nexttile(layout);
    plot(ax, horizons, modeIRF(:, series), "k", LineWidth = 1.2);
    hold(ax, "on");
    plot(ax, horizons, squeeze(bandIRF(:, series, :)), "-.", Color = [0.75 0.1 0.1]);
    yline(ax, 0, ":", Color = [0.45 0.45 0.45]);
    xlim(ax, [horizons(1) horizons(end)]);
    title(ax, seriesNames(series)); box(ax, "on");
end

legend(ax, ["IRF at mode", compose("%gth quantile", quantiles*100)], ...
    Location = "northeastoutside");
title(layout, plotTitle);
end

function relative = localRelativeNorm(x, y)
relative = norm(x(:) - y(:))/max(norm(y(:)), eps);
end
