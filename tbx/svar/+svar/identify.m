function [impact, info] = identify(Mdl, specification, nvp)
%IDENTIFY - Identify structural shocks from a reduced-form covariance
%   IMPACT = IDENTIFY(Mdl,METHOD) identifies structural shocks using the
%   covariance matrix in the fitted VARM model Mdl. For Bayesian VAR models,
%   specify the covariance draw using Sigma.
%
%   [IMPACT,INFO] = IDENTIFY(...) also returns the resolved identification
%   method, ordering, shock names, and number of covariance draws.
%
%   IMPACT = IDENTIFY(...,Sigma=SIGMA) identifies one covariance matrix or
%   an m-by-m-by-N stack of posterior covariance draws. The returned IMPACT
%   has the matching m-by-m or m-by-m-by-N size.
%
%   IMPACT = IDENTIFY(...,Ordering=ORDERING) uses the supplied complete
%   recursive ordering. ORDERING is either a permutation of series indices
%   or of the model's series names. It determines the structural-shock
%   column order; response rows remain in the model's native series order.
%
%   IMPACT = IDENTIFY(Mdl,RESTRICTIONS,...) draws structural models that
%   satisfy public svar restriction declarations. This form requires
%   ShockNames and uses independent reduced-form draws. For impact zero
%   restrictions, rotations are drawn conditionally on the zero constraints.
%   Set Sampler="elliptical" to use an elliptical slice within Gibbs sampler
%   for sign-restricted conjugate BVAR posteriors.
%
%   METHOD must be "cholesky", "recursive", "triangular", or "reduced".
%   Cholesky and recursive identification return a unit-variance impact
%   matrix. Triangular identification returns the same impact matrix and
%   additionally reports BEAR's unit-diagonal structural matrix and diagonal
%   structural-shock covariance in INFO. Reduced identification returns an
%   identity impact matrix without orthogonalizing reduced-form innovations.
%
%   See also svar.irf, svar.fevd, historicalDecomposition, chol

arguments
    Mdl (1,1) {mustBeA(Mdl, ["varm", "bvar"])}
    specification = "cholesky"
    nvp.Sigma double {mustBeReal, mustBeFinite} = []
    nvp.Ordering = []
    nvp.ShockNames (1,:) string = strings(1, 0)
    nvp.NumDraws (1,1) double {mustBeInteger, mustBePositive} = 1
    nvp.BatchSize (1,1) double {mustBeInteger, mustBePositive} = 100
    nvp.MaxCandidates (1,1) double {mustBeInteger, mustBePositive} = 100000
    nvp.Sampler (1,1) string {mustBeMember(nvp.Sampler, ...
        ["independent", "elliptical"])} = "independent"
    nvp.BurnIn (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    nvp.Thin (1,1) double {mustBeInteger, mustBePositive} = 1
    nvp.MaxSliceSteps (1,1) double {mustBeInteger, mustBePositive} = 10000
    nvp.NumParticles (1,1) double {mustBeInteger, mustBePositive} = 500
    nvp.MutationSteps (1,1) double {mustBeInteger, mustBePositive} = 2
    nvp.TargetSurvivalRate (1,1) double { ...
        mustBeGreaterThan(nvp.TargetSurvivalRate, 0), ...
        mustBeLessThanOrEqual(nvp.TargetSurvivalRate, 1)} = 0.2
    nvp.MaxInitializationStages (1,1) double {mustBeInteger, ...
        mustBePositive} = 50
end

if isa(specification, "svar.identification.Restriction")
    if isempty(nvp.ShockNames)
        error("svar:identify:MissingShockNames", ...
            "Restricted identification requires complete ShockNames.");
    end
    if numel(nvp.ShockNames) ~= Mdl.NumSeries || ...
            numel(unique(nvp.ShockNames)) ~= Mdl.NumSeries
        error("svar:identify:InvalidShockNames", ...
            "ShockNames must contain one unique name for each series.");
    end
    restrictionSet = svar.identification.RestrictionSet( ...
        string(Mdl.SeriesNames), ShockOrder=nvp.ShockNames);
    for restriction = reshape(specification, 1, [])
        restrictionSet = restrictionSet.add(restriction);
    end
    switch nvp.Sampler
        case "independent"
            [impact, info] = svar.identification.IndependentSampler.identify( ...
                Mdl, restrictionSet, ...
                NumDraws=nvp.NumDraws, BatchSize=nvp.BatchSize, ...
                MaxCandidates=nvp.MaxCandidates);
        case "elliptical"
            if ~isa(Mdl, "conjugatebvarm")
                error("svar:identify:EllipticalRequiresConjugate", ...
                    "Elliptical sampling requires a conjugatebvarm posterior.");
            end
            [impact, info] = svar.identification.EllipticalSampler.identify( ...
                Mdl, restrictionSet, NumDraws=nvp.NumDraws, ...
                BurnIn=nvp.BurnIn, Thin=nvp.Thin, ...
                MaxSliceSteps=nvp.MaxSliceSteps, ...
                NumParticles=nvp.NumParticles, ...
                MutationSteps=nvp.MutationSteps, ...
                TargetSurvivalRate=nvp.TargetSurvivalRate, ...
                MaxInitializationStages=nvp.MaxInitializationStages);
    end
    return
end

if ~isstring(specification) || ~isscalar(specification)
    error("svar:identify:InvalidSpecification", ...
        "Specify a method string or svar restriction declarations.");
end

switch specification
    case {"cholesky", "recursive"}
        method = "cholesky";
    case "triangular"
        method = "triangular";
    case "reduced"
        method = "reduced";
    otherwise
        error("svar:identify:UnsupportedMethod", ...
            "Method must be ""cholesky"", ""recursive"", ""triangular"", or ""reduced"".");
end
        
% Grab series
numSeries = Mdl.NumSeries;
seriesNames = Mdl.SeriesNames;

% Check variable ordering
ordering = localResolveOrdering(nvp.Ordering, seriesNames);

% Select covariance
if ~isempty(nvp.Sigma)
    sigma = nvp.Sigma;
elseif isa(Mdl, "varm")
    sigma = Mdl.Covariance;
else
    error("svar:identify:MissingCovarianceDraw", ...
        "For a Bayesian VAR, specify Sigma as a covariance draw or stack of draws.");
end
numDraws = localValidateCovariance(sigma, numSeries);

% Initialize Impact Matrix
if ismatrix(sigma)
    impact = zeros(numSeries, numSeries);
    structuralMatrix = zeros(numSeries, numSeries);
    shockCovariance = zeros(numSeries, numSeries);
else
    impact = zeros(numSeries, numSeries, numDraws);
    structuralMatrix = zeros(numSeries, numSeries, numDraws);
    shockCovariance = zeros(numSeries, numSeries, numDraws);
end

sigmaOrdered = sigma(ordering, ordering, :);
for draw = 1:numDraws
    covariance = sigmaOrdered(:, :, draw);
    switch method
        case "cholesky"
            [factor, matrix, shockCov] = ...
                svar.identification.choleskyIdentification(covariance, draw);
        case "triangular"
            [factor, matrix, shockCov] = ...
                svar.identification.triangularIdentification(covariance, draw);
        case "reduced"
            [factor, matrix, shockCov] = ...
                svar.identification.reducedIdentification(covariance);
    end
    impact(ordering, :, draw) = factor;
    structuralMatrix(ordering, :, draw) = matrix;
    shockCovariance(:, :, draw) = shockCov;
end

info = struct( ...
    "Method", method, ...
    "Ordering", ordering, ...
    "ShockOrder", seriesNames(ordering), ...
    "NumDraws", numDraws, ...
    "StructuralMatrix", structuralMatrix, ...
    "ShockCovariance", shockCovariance);
end

function ordering = localResolveOrdering(ordering, seriesNames)
numSeries = numel(seriesNames);
if isempty(ordering)
    ordering = 1:numSeries;
elseif isnumeric(ordering)
    if ~isvector(ordering) || ~isreal(ordering) || any(~isfinite(ordering)) ...
            || any(ordering ~= fix(ordering))
        localInvalidOrdering();
    end
    ordering = reshape(ordering, 1, []);
    if ~isequal(sort(ordering), 1:numSeries)
        localInvalidOrdering();
    end
elseif isstring(ordering) || iscellstr(ordering)
    orderNames = reshape(string(ordering), 1, []);
    [isMember, ordering] = ismember(orderNames, seriesNames);
    if numel(ordering) ~= numSeries || ~all(isMember) ...
            || numel(unique(ordering)) ~= numSeries
        localInvalidOrdering();
    end
else
    localInvalidOrdering();
end
end

function localInvalidOrdering()
error("svar:identify:InvalidOrdering", ...
    "Ordering must be a complete, unique permutation of series names or indices.");
end

function numDraws = localValidateCovariance(sigma, numSeries)
if ~isnumeric(sigma) || ~isreal(sigma) || any(~isfinite(sigma), "all") ...
        || ~(ismatrix(sigma) || ndims(sigma) == 3) ...
        || size(sigma, 1) ~= numSeries || size(sigma, 2) ~= numSeries
    error("svar:identify:InvalidCovariance", ...
        "Sigma must be a finite %d-by-%d matrix or a stack of such matrices.", ...
        numSeries, numSeries);
end
numDraws = size(sigma, 3);
end
