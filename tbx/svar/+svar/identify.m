function [impact, info] = identify(Mdl, Method, nvp)
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
    Method (1,1) string = "cholesky"
    nvp.Sigma double {mustBeReal, mustBeFinite} = []
    nvp.Ordering = []
end

switch Method
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
