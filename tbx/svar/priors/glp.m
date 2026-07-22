function [mdl, info] = glp(numseries, numlags, Y, nvp, nvp2)
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

    function negObj = objective(x)
        [lambdaArgs, candPsi] = unpackAll(x);
        try
            mdlCandidate = minnesotamniwbvarm(numseries, numlags, Y, ...
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
mdl = minnesotamniwbvarm(numseries, numlags, Y, ...
    "Psi", finalPsi, finalLambdaArgs{:}, fixedArgs{:}, buildArgs{:});

if nargout > 1
    info = struct( ...
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
