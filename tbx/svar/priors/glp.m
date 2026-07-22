function [mdl, info, chain] = glp(numseries, numlags, Y, nvp, nvp2)
%GLP Tune Minnesota hyperparameters by marginal likelihood.
%   mdl = GLP(NUMSERIES,NUMLAGS,Y) tunes the free hyperparameters of a
%   MINNESOTAMNIWBVARM prior by maximising its analytic marginal likelihood
%   against Y, and returns the prior built at the maximiser.
%
%   Each hyperparameter is given in one of three states:
%       scalar        FIXED at that value, not tuned
%       [lower upper] FREE, tuned within those bounds under a flat prior
%       hyperprior    FREE, tuned within the object's Bounds under its log prior
%   So the free/fixed split is implicit in what you pass - there is no
%   separate list of free parameters to keep in sync.
%
%   mdl = GLP(...,lambda1=L1,lambda3=L3,lambda4=L4,lambda5=L5,Vc=VC,
%   PriorMean=PM) sets the Minnesota hyperparameters. lambda4 and lambda5
%   default to Inf, which switches the sum-of-coefficients and dummy-initial-
%   observation priors off. Vc and PriorMean are always fixed.
%
%   mdl = GLP(...,Psi=PSI) controls the residual-variance scale:
%       "exact" | "conditional"   estimated from Y once, then held fixed
%       1-by-n numeric            fixed at these variances
%       2-by-n numeric            free within these [lower; upper] bounds
%       scalar hyperprior         free, broadcast independently to all n
%       1-by-n hyperprior         free, one log prior per series
%   The default is "exact". A string Psi is resolved ONCE before the search,
%   never per candidate.
%
%   [mdl,info] = GLP(...) also returns a struct describing the search: which
%   hyperparameters were free, their bounds and starting point, the
%   maximiser, and the fmincon exit information.
%
%   Integrating over the hyperparameters
%   ------------------------------------
%   [mdl,info,chain] = GLP(...,NumDraws=M) additionally samples the
%   hyperparameter posterior instead of only maximising it, and returns M
%   draws of the hyperparameters AND of the VAR parameters they imply. This
%   is the Metropolis step of Giannone, Lenza and Primiceri (2012), appendix
%   B: the maximiser above is only its step 1.
%
%   The sampler is a random walk whose proposal covariance is the inverse
%   Hessian of the negative log posterior at the maximiser, computed by
%   central differences, scaled by ProposalScale squared. Because the search
%   runs in natural, box-constrained coordinates, that Hessian needs no
%   reparameterisation; a candidate outside the box is simply rejected. Tune
%   ProposalScale for an acceptance rate around 0.2-0.3 and check
%   chain.AcceptanceRate.
%
%       NumDraws       draws KEPT (0, the default, means maximise only)
%       BurnIn         draws discarded first (default: NumDraws)
%       ProposalScale  random-walk step size c (default 1)
%
%   Bands built from chain widen those built from mdl alone, because the
%   latter condition on one value of the hyperparameters.
%
%   Example
%   -------
%       mdl = glp(size(Y,2), 4, Y, ...
%           lambda1 = hyperprior("Gamma",0.2,0.4,Bounds=[1e-4 5]), ...
%           lambda4 = hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
%           lambda5 = hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
%           Psi     = "conditional");
%
%   See also MINNESOTAMNIWBVARM, HYPERPRIOR, LOGMARGINALLIKELIHOOD.

arguments
    numseries (1,1) double {mustBeInteger, mustBePositive}
    numlags   (1,1) double {mustBeInteger, mustBePositive}
    Y               {mustBeNonempty}
    nvp.lambda1     = 0.2
    nvp.lambda3     = 1
    nvp.lambda4     = Inf
    nvp.lambda5     = Inf
    nvp.Vc        (1,1) double {mustBePositive} = 1e4
    nvp.PriorMean (1,:) double = 1
    nvp.Psi         = "exact"
    nvp.OptimOptions      = optimoptions("fmincon", ...
        Display = "final", ...
        FiniteDifferenceStepSize = 1e-4, ...
        FunctionTolerance = 1e-12, ...
        StepTolerance = 1e-12, ...
        ConstraintTolerance = 1e-12)
    nvp.NumDraws      (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    nvp.BurnIn              double {mustBeScalarOrEmpty, mustBeInteger, ...
                                    mustBeNonnegative} = []
    nvp.ProposalScale (1,1) double {mustBePositive} = 1
    nvp2.Description
    nvp2.IncludeConstant
    nvp2.IncludeTrend
    nvp2.NumPredictors
    nvp2.SeriesNames
end

if istabular(Y)
    if ~isfield(nvp2, "SeriesNames")
        nvp2.SeriesNames = string(Y.Properties.VariableNames);
    end
    Y = Y{:,:};
end

if ~isnumeric(Y)
    error("glp:BadDataType", "Responses input must be numeric or tabular.")
end

if size(Y, 2) ~= numseries
    error("glp:invalidData", ...
        "Y must have %d columns, one per series.", numseries);
end

buildArgs = namedargs2cell(nvp2);
fixedArgs = {"Vc", nvp.Vc, "PriorMean", nvp.PriorMean};

lambdaNames = ["lambda1","lambda3","lambda4","lambda5"];
[lamX0, lamLB, lamUB, lamFixed, lamPrior, freeLambdas] = ...
    localPackLambdas(nvp, lambdaNames);
[psiX0, psiLB, psiUB, psiFixed, psiPrior, psiNames] = ...
    localPackPsi(nvp.Psi, numseries, numlags, Y);

x0 = [lamX0, psiX0];
lb = [lamLB, psiLB];
ub = [lamUB, psiUB];
nLambdaFree = numel(freeLambdas);
psiFree = ~isempty(psiX0);
useHyperprior = ~isempty(fieldnames(lamPrior)) || ~isempty(psiPrior);

    function [lambdaArgs, candPsi] = unpackAll(x)
        lambdaValues = lamFixed;
        for k = 1:nLambdaFree
            lambdaValues.(freeLambdas(k)) = x(k);
        end
        lambdaArgs = namedargs2cell(lambdaValues);

        if psiFree
            candPsi = x(nLambdaFree+1:end);
        else
            candPsi = psiFixed;
        end
    end

    function [negObj, mdlCandidate] = objective(x)
        % The candidate prior is returned as well so the sampler can draw
        % (Coeff, Sigma) from it without rebuilding what was just built.
        [lambdaArgs, candPsi] = unpackAll(x);
        mdlCandidate = [];
        try
            mdlCandidate = svar.minnesotamniwbvarm(numseries, numlags, Y, ...
                "Psi", candPsi, lambdaArgs{:}, fixedArgs{:}, buildArgs{:});
            negObj = -logMarginalLikelihood(mdlCandidate);
        catch
            negObj = 1e10;
        end

        if isfinite(negObj) && negObj <= 1e10
            logPrior = localLogHyperprior(lamPrior, x(1:nLambdaFree), freeLambdas) ...
                + localLogPsiHyperprior(psiPrior, candPsi);
            negObj = negObj - logPrior;
        end

        if ~isfinite(negObj) || negObj > 1e10
            negObj = 1e10;
        end
    end

if isempty(x0)
    xHat = x0;
    fminconOutput = struct("Skipped", true, ...
        "Reason", "No free hyperparameters.");
    fval = objective(xHat);
    exitflag = NaN;
else
    [xHat, fval, exitflag, fminconOutput] = fmincon(@objective, x0, ...
        [], [], [], [], lb, ub, [], nvp.OptimOptions);
end

[finalLambdaArgs, finalPsi] = unpackAll(xHat);
mdl = svar.minnesotamniwbvarm(numseries, numlags, Y, ...
    "Psi", finalPsi, finalLambdaArgs{:}, fixedArgs{:}, buildArgs{:});

if nvp.NumDraws > 0
    if isempty(x0)
        error("glp:nothingToSample", ...
            "NumDraws requires at least one free hyperparameter.");
    end
    burnIn = nvp.BurnIn;
    if isempty(burnIn)
        burnIn = nvp.NumDraws;
    end
    hessian = localNumericalHessian(@objective, xHat, lb, ub);
    chain = localMetropolis(@objective, @unpackAll, xHat, -fval, hessian, ...
        lb, ub, nvp.NumDraws, burnIn, nvp.ProposalScale, ...
        lamFixed, freeLambdas, lambdaNames, numseries);
else
    hessian = [];
    chain   = struct([]);
end

if nargout > 1
    info = struct( ...
        "NumDraws",        nvp.NumDraws, ...
        "ProposalScale",   nvp.ProposalScale, ...
        "Hessian",         hessian, ...
        "FreeLambdas",     freeLambdas, ...
        "PsiNames",        psiNames, ...
        "PsiFree",         psiFree, ...
        "FinalLambdas",    localLambdaStruct(mdl, lambdaNames), ...
        "InitialPsi",      nvp.Psi, ...
        "FinalPsi",        finalPsi, ...
        "UsedHyperprior",  useHyperprior, ...
        "X0",              x0, ...
        "LowerBound",      lb, ...
        "UpperBound",      ub, ...
        "XHat",            xHat, ...
        "Objective",       fval, ...
        "ExitFlag",        exitflag, ...
        "FminconOutput",   fminconOutput);
end

end

function chain = localMetropolis(objective, unpackAll, xHat, logPostHat, ...
    hessian, lb, ub, numDraws, burnIn, scale, lamFixed, freeLambdas, ...
    lambdaNames, numseries)
%LOCALMETROPOLIS Random-walk Metropolis over the free hyperparameters.
%   Giannone, Lenza and Primiceri (2012), appendix B. The chain starts AT the
%   maximiser (their step 1), proposes from a Gaussian centred on the current
%   draw with covariance scale^2 * inv(hessian), and accepts on the log
%   posterior ratio. The proposal is symmetric, so no Hastings correction.
%
%   Conditional on each accepted hyperparameter draw, (Coeff, Sigma) come
%   from the exact Normal-Inverse-Wishart posterior (their step 4) - no inner
%   chain, so the kept draws are independent given the hyperparameters.

proposalCovariance = localProposalCovariance(hessian, scale);

totalDraws = burnIn + numDraws;
x          = xHat;
logPost    = logPostHat;

% The posterior for the CURRENT hyperparameters is cached: a rejected step
% leaves the hyperparameters unchanged, so it would rebuild the same object.
[~, currentPrior] = objective(x);
currentPosterior  = estimate(currentPrior, Display = "off");
numCoefficients   = numel(currentPosterior.Mu)/numseries;

coefficientDraws = zeros(numCoefficients, numseries, numDraws);
covarianceDraws  = zeros(numseries, numseries, numDraws);
xDraws           = zeros(numDraws, numel(xHat));
logPostDraws     = zeros(numDraws, 1);
accepted         = 0;

for draw = 1:totalDraws
    candidate = mvnrnd(x, proposalCovariance);

    % Box constraints are enforced by rejection: fmincon searched a box, and
    % the posterior is only defined inside it.
    if all(candidate >= lb) && all(candidate <= ub)
        [negObj, candidatePrior] = objective(candidate);
        candidateLogPost = -negObj;

        if log(rand) < candidateLogPost - logPost
            x                = candidate;
            logPost          = candidateLogPost;
            currentPosterior = estimate(candidatePrior, Display = "off");
            accepted         = accepted + 1;
        end
    end

    if draw > burnIn
        kept = draw - burnIn;
        [coefficients, sigma] = simulate(currentPosterior, NumDraws = 1);

        coefficientDraws(:,:,kept) = reshape(coefficients, [], numseries);
        covarianceDraws(:,:,kept)  = sigma;
        xDraws(kept,:)             = x;
        logPostDraws(kept)         = logPost;
    end
end

chain = struct( ...
    "Coefficients",       coefficientDraws, ...
    "Sigma",              covarianceDraws, ...
    "LogPosterior",       logPostDraws, ...
    "AcceptanceRate",     accepted/totalDraws, ...
    "NumDraws",           numDraws, ...
    "BurnIn",             burnIn, ...
    "ProposalScale",      scale, ...
    "ProposalCovariance", proposalCovariance, ...
    "FreeLambdas",        freeLambdas, ...
    "X",                  xDraws);

% Expand the free vector back into named hyperparameters, so a consumer never
% has to know the packing order.
psiDraws = zeros(numDraws, numseries);
for name = lambdaNames
    chain.(name) = zeros(numDraws, 1);
end
for kept = 1:numDraws
    [lambdaArgs, psi] = unpackAll(xDraws(kept,:));
    lambdaValues      = localArgsToStruct(lambdaArgs, lamFixed);
    for name = lambdaNames
        chain.(name)(kept) = lambdaValues.(name);
    end
    psiDraws(kept,:) = psi;
end
chain.Psi = psiDraws;
end

function hessian = localNumericalHessian(objective, x, lb, ub)
%LOCALNUMERICALHESSIAN Central-difference Hessian of the negative log posterior.
%   Deliberately NOT fmincon's 7th output: that is a quasi-Newton
%   approximation to the Hessian of the Lagrangian, contaminated by the
%   barrier terms, and on this objective it comes back with a condition
%   number around 1e12 - inverting it produces proposals hundreds of times
%   wider than the parameters themselves and an acceptance rate of zero.
%
%   Step size matters more than usual because second differences divide by
%   h^2: a step tuned for gradients is swamped by roundoff here. Around 1% of
%   each parameter sits on the plateau where truncation and roundoff are both
%   small; below ~0.1% the estimate diverges.
n    = numel(x);
step = 1e-2*abs(x);
step(step == 0) = 1e-2;

% Every evaluation point must stay inside the box the search used.
step = min(step, 0.9*min(x - lb, ub - x));
if any(step <= 0)
    error("glp:degenerateHessian", ...
        "Cannot take a finite-difference step inside the bounds; a " + ...
        "hyperparameter sits on its bound.");
end

    function value = at(shift)
        value = objective(x + shift);
    end

f0      = at(zeros(size(x)));
hessian = zeros(n);
e       = @(k) double((1:n) == k);

forward  = arrayfun(@(k) at( step(k)*e(k)), 1:n);
backward = arrayfun(@(k) at(-step(k)*e(k)), 1:n);

for a = 1:n
    hessian(a,a) = (forward(a) - 2*f0 + backward(a))/step(a)^2;
    for b = a+1:n
        shiftA = step(a)*e(a);
        shiftB = step(b)*e(b);
        cross  = at(shiftA + shiftB) - at(shiftA - shiftB) ...
               - at(-shiftA + shiftB) + at(-shiftA - shiftB);
        hessian(a,b) = cross/(4*step(a)*step(b));
        hessian(b,a) = hessian(a,b);
    end
end

hessian = (hessian + hessian')/2;
end

function proposalCovariance = localProposalCovariance(hessian, scale)
%LOCALPROPOSALCOVARIANCE Inverse Hessian, forced symmetric positive definite.
%   A finite-difference Hessian can still be slightly indefinite in flat
%   directions. Reflecting the eigenvalues (as the authors' own code does)
%   keeps the proposal usable instead of failing at the first draw.
hessian = (hessian + hessian')/2;

[vectors, values] = eig(hessian);
values            = abs(diag(values));
values(values < eps(max(values))) = eps(max(values));

precision          = vectors*diag(values)*vectors';
proposalCovariance = (scale^2)*(precision\eye(size(precision)));
proposalCovariance = (proposalCovariance + proposalCovariance')/2;
end

function lambdaValues = localArgsToStruct(lambdaArgs, lamFixed)
lambdaValues = lamFixed;
for i = 1:2:numel(lambdaArgs)
    lambdaValues.(lambdaArgs{i}) = lambdaArgs{i+1};
end
end

function [x0, lb, ub, fixedValues, priors, freeNames] = localPackLambdas(nvp, names)
%LOCALPACKLAMBDAS Split the lambdas into fixed values and a free search box.
x0 = [];
lb = [];
ub = [];
freeNames = string.empty(1,0);
priors = struct();
fixedValues = struct();

for i = 1:numel(names)
    name = names(i);
    value = nvp.(name);

    if isa(value, "hyperprior")
        if ~isscalar(value)
            error("glp:invalidHyperparameter", ...
                "%s must be a scalar hyperprior.", name);
        end
        freeNames(end+1) = name; %#ok<AGROW>
        priors.(name) = value;
        x0(end+1) = value.X0; %#ok<AGROW>
        lb(end+1) = value.Bounds(1); %#ok<AGROW>
        ub(end+1) = value.Bounds(2); %#ok<AGROW>
        continue
    end

    if ~isnumeric(value) || ~isvector(value) || numel(value) > 2
        error("glp:invalidHyperparameter", ...
            "%s must be a scalar, a 2-element [lower upper] bound, " + ...
            "or a hyperprior.", name);
    end

    if any(value <= 0)
        error("glp:invalidHyperparameter", ...
            "%s must be positive.", name);
    end

    if isscalar(value)
        fixedValues.(name) = value;
        continue
    end

    if any(~isfinite(value)) || value(1) >= value(2)
        error("glp:invalidHyperparameter", ...
            "%s bounds must be finite and satisfy lower < upper.", name);
    end
    freeNames(end+1) = name; %#ok<AGROW>
    x0(end+1) = sqrt(value(1)*value(2)); %#ok<AGROW>
    lb(end+1) = value(1); %#ok<AGROW>
    ub(end+1) = value(2); %#ok<AGROW>
end
end

function [x0, lb, ub, fixedPsi, psiPrior, names] = localPackPsi(Psi, numseries, numlags, Y)
names = "Psi" + string(1:numseries);
fixedPsi = [];
psiPrior = [];

if isstring(Psi) || ischar(Psi)
    % Resolve the estimator ONCE: the string form refits an AR per series,
    % which must never happen inside the objective.
    fixedPsi = estimateResidualVariances(Y, numlags, Method=string(Psi));
    x0 = [];
    lb = [];
    ub = [];
    return
end

if isa(Psi, "hyperprior")
    psiPrior = reshape(Psi, 1, []);
    if isscalar(psiPrior) && numseries > 1
        psiPrior = repmat(psiPrior, 1, numseries);
    end
    if numel(psiPrior) ~= numseries
        error("glp:invalidPsi", ...
            "Psi hyperprior input must be scalar or have one element per series.");
    end
    x0 = arrayfun(@(p) p.X0, psiPrior);
    lb = arrayfun(@(p) p.Bounds(1), psiPrior);
    ub = arrayfun(@(p) p.Bounds(2), psiPrior);
    return
end

if ~isnumeric(Psi)
    error("glp:invalidPsi", ...
        "Psi must be an estimator name, numeric, or a hyperprior array.");
end

if isvector(Psi) && numel(Psi) == numseries
    fixedPsi = reshape(Psi, 1, []);
    if any(~isfinite(fixedPsi)) || any(fixedPsi <= 0)
        error("glp:invalidPsi", ...
            "Fixed Psi values must be finite and positive.");
    end
    x0 = [];
    lb = [];
    ub = [];
    return
end

if ismatrix(Psi) && isequal(size(Psi), [2 numseries])
    lb = Psi(1,:);
    ub = Psi(2,:);
    if any(~isfinite(lb)) || any(~isfinite(ub)) || any(lb <= 0) || any(lb >= ub)
        error("glp:invalidPsi", ...
            "Psi bounds must be finite, positive, and ordered lower < upper.");
    end
    x0 = sqrt(lb .* ub);
    return
end

error("glp:invalidPsi", ...
    "Psi must be an estimator name, 1-by-n numeric, 2-by-n numeric bounds, " + ...
    "scalar hyperprior, or 1-by-n hyperprior.");
end

function lp = localLogHyperprior(priors, x, names)
lp = 0;
for i = 1:numel(names)
    if isfield(priors, names(i))
        lp = lp + priors.(names(i)).logpdf(x(i));
    end
end
end

function lp = localLogPsiHyperprior(psiPrior, psi)
lp = 0;
if isempty(psiPrior)
    return
end
for i = 1:numel(psi)
    lp = lp + psiPrior(i).logpdf(psi(i));
end
end

function s = localLambdaStruct(mdl, names)
s = struct();
for i = 1:numel(names)
    s.(names(i)) = mdl.(names(i));
end
end
