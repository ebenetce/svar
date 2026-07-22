function HD = historicalDecomposition(Mdl, Impact, Y, opts)
%historicalDecomposition Historical decomposition for structural VAR models.
%   HD = historicalDecomposition(Mdl, Impact, Y) decomposes the realized
%   history in Y into the cumulative contributions of each structural shock.
%   Mdl can be either a fully specified VARM model or a Bayesian VAR model.
%   When Mdl is a Bayesian VAR object, the function converts it to an
%   equivalent VARM model using BVAR2VAR before inferring reduced-form
%   residuals.
%
%   HD = historicalDecomposition(Mdl, Impact, Y, Name=Value) controls the
%   decomposition:
%       Y0             - Presample responses for infer
%       X              - Predictor data for infer
%       Residuals      - Precomputed reduced-form residuals
%       VariableIndices - Variables to report
%       ShockIndices    - Shocks to report
%
%   Output HD is a struct with fields:
%       Model                - VARM model used for the decomposition
%       Impact               - Structural impact matrix
%       Residuals            - Reduced-form residuals
%       StructuralShocks     - Structural shocks corresponding to ShockIndices
%       FullContributions    - T-by-K-by-K contributions (time, variable, shock)
%       Contributions        - T-by-numVariables-by-numShocks selected subset
%       Total                - T-by-numVariables total contribution of selected shocks
%       VariableIndices      - Selected variable indices
%       ShockIndices         - Selected shock indices
%       VariableNames        - Selected variable names
%       ShockNames           - Selected shock names
%
%   Example
%       EstMdl = estimate(varm(2,1), Y, Display="off");
%       Impact = chol(EstMdl.Covariance, "lower");
%       HD = historicalDecomposition(EstMdl, Impact, Y);
%
%   See also infer, varm, bvar2var

arguments
    Mdl (1,1) %{mustBeVarOrBvar}
    Impact (:,:) double
    Y = []
    opts.Y0 double = []
    opts.X double = []
    opts.Residuals double = []
    opts.VariableIndices = []
    opts.ShockIndices = []
end

varMdl = bvar2var(Mdl);
numSeries = Mdl.NumSeries;

if size(Impact, 1) ~= numSeries || size(Impact, 2) ~= numSeries
    error("historicalDecomposition:InvalidImpactSize", ...
        "Impact must be a %d-by-%d matrix.", numSeries, numSeries);
end

if isempty(opts.Residuals)
    if isempty(Y)
        error("historicalDecomposition:MissingData", ...
            "Y is required unless Residuals is provided.");
    end

    nvArgs = {};
    if ~isempty(opts.Y0), nvArgs = [nvArgs, {'Y0', opts.Y0}]; end
    if ~isempty(opts.X),  nvArgs = [nvArgs, {'X', opts.X}];  end
    residuals = infer(varMdl, Y, nvArgs{:});
else
    residuals = opts.Residuals;
end

if size(residuals, 2) ~= numSeries
    error("historicalDecomposition:InvalidResidualSize", ...
        "Residuals must have %d columns.", numSeries);
end

variableIndices = localResolveIndices(opts.VariableIndices, numSeries);
shockIndices = localResolveIndices(opts.ShockIndices, numSeries);

numObs = size(residuals, 1);
irf = localImpulseResponses(varMdl, Impact, numObs - 1);
structuralShocks = (Impact \ residuals.').';
fullContributions = localContributions(irf, structuralShocks);

seriesNames = string(varMdl.SeriesNames);
if isempty(seriesNames)
    seriesNames = "Series" + (1:numSeries);
end

HD = struct();
HD.Model = varMdl;
HD.Impact = Impact;
HD.Residuals = residuals;
HD.StructuralShocks = structuralShocks(:, shockIndices);
HD.FullContributions = fullContributions;
HD.Contributions = fullContributions(:, variableIndices, shockIndices);
HD.Total = sum(HD.Contributions, 3);
HD.VariableIndices = variableIndices;
HD.ShockIndices = shockIndices;
HD.VariableNames = seriesNames(variableIndices);
HD.ShockNames = seriesNames(shockIndices);
end

function indices = localResolveIndices(indices, upperBound)
if isempty(indices)
    indices = 1:upperBound;
else
    mustBePositive(indices);
    mustBeInteger(indices);
    if any(indices > upperBound)
        error("historicalDecomposition:InvalidIndex", ...
            "Indices must be between 1 and %d.", upperBound);
    end
    indices = reshape(indices, 1, []);
end
end

function irf = localImpulseResponses(varMdl, impact, horizon)
% AR: 1-by-p cell array of numSeries-by-numSeries AR coefficient matrices (Mdl.AR)
numSeries = size(impact,1);
p = varMdl.P;
AR = varMdl.AR;
Theta = cell(horizon+1,1);
Theta{1} = eye(numSeries);
for h = 1:horizon
    acc = zeros(numSeries);
    for j = 1:min(h,p)
        acc = acc + AR{j} * Theta{h-j+1};
    end
    Theta{h+1} = acc;
end
irf = zeros(numSeries^2, horizon+1);
for h = 0:horizon
    theta = Theta{h+1} * impact;
    irf(:, h+1) = theta(:);
end
end

function contributions = localContributions(irf, structuralShocks)
numObs = size(structuralShocks, 1);
numSeries = size(structuralShocks, 2);
contributions = zeros(numObs, numSeries, numSeries);

for shock = 1:numSeries
    rows = (shock - 1) * numSeries + (1:numSeries);
    x = structuralShocks(:, shock);
    for variable = 1:numSeries
        b = irf(rows(variable), 1:numObs);
        contributions(:, variable, shock) = filter(b, 1, x);
    end
end
end
