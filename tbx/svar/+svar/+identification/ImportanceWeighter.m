classdef ImportanceWeighter
    %IMPORTANCEWEIGHTER - ARW importance weights for impact zero restrictions

    properties (SetAccess = private)
        Info (1,1) struct
    end

    methods
        function obj = ImportanceWeighter(restrictions, numCoefficients)
            obj.Info = localBuildInfo(restrictions, numCoefficients);
        end

        function weight = logWeight(obj, coefficients, sigma, rotation)
            info = obj.Info;
            parameters = localFHInv([coefficients(:); sigma(:); rotation(:)], ...
                info);
            A0 = reshape(parameters(1:info.NumSeries^2), ...
                info.NumSeries, info.NumSeries);
            structuralJacobian = 0.5*info.NumSeries*(info.NumSeries + 1)*log(2) ...
                - (2*info.NumSeries + info.NumCoefficients + 1)* ...
                localLogAbsDet(A0);
            volume = localLogVolume(@(x) localFFH(x, info), parameters, ...
                @(x) localZeroRestrictions(x, info));
            weight = structuralJacobian - volume;
        end
    end
end

function info = localBuildInfo(restrictions, numCoefficients)
numSeries = numel(restrictions.SeriesNames);
selectors = cell(numSeries, 1);
weights = cell(numSeries, 1);
dimension = 0;
numZeros = 0;
for shock = 1:numSeries
    selectors{shock} = restrictions.zeroSelectorMatrix(shock);
    degreesOfFreedom = numSeries - (shock - 1 + size(selectors{shock}, 1));
    if degreesOfFreedom < 0
        error("svar:identify:OverIdentifiedZeroRestrictions", ...
            "Shock ""%s"" has too many impact zero restrictions.", ...
            restrictions.ShockOrder(shock));
    end
    weights{shock} = randn(degreesOfFreedom, numSeries);
    dimension = dimension + degreesOfFreedom;
    numZeros = numZeros + size(selectors{shock}, 1);
end
info = struct( ...
    "NumSeries", numSeries, ...
    "NumCoefficients", numCoefficients, ...
    "Selectors", {selectors}, ...
    "Weights", {weights}, ...
    "Dimension", dimension, ...
    "NumZeros", numZeros);
end

function y = localFHInv(x, info)
numSeries = info.NumSeries;
numCoefficients = info.NumCoefficients;
B = reshape(x(1:numCoefficients*numSeries), numCoefficients, numSeries);
sigma = reshape(x(numCoefficients*numSeries + 1: ...
    (numCoefficients + numSeries)*numSeries), numSeries, numSeries);
rotation = reshape(x((numCoefficients + numSeries)*numSeries + 1:end), ...
    numSeries, numSeries);
A0 = chol((sigma + sigma.')/2) \ rotation;
Aplus = B*A0;
y = [A0(:); Aplus(:)];
end

function y = localFH(x, info)
numSeries = info.NumSeries;
numCoefficients = info.NumCoefficients;
A0 = reshape(x(1:numSeries^2), numSeries, numSeries);
Aplus = reshape(x(numSeries^2 + 1:end), numCoefficients, numSeries);
B = Aplus/A0;
sigma = inv(A0*A0.');
rotation = chol((sigma + sigma.')/2)*A0;
y = [B(:); sigma(:); rotation(:)];
end

function y = localFFH(x, info)
numSeries = info.NumSeries;
z = localFH(x, info);
identity = eye(numSeries);
zeroResponses = localImpactResponses(localFHInv( ...
    [z(1:numSeries*(info.NumCoefficients + numSeries)); identity(:)], info), ...
    info);
rotation = reshape(z(numSeries*(info.NumCoefficients + numSeries) + 1:end), ...
    numSeries, numSeries);
w = zeros(info.Dimension, 1);
offset = 0;
for shock = 1:numSeries
    count = size(info.Weights{shock}, 1);
    matrix = [rotation(:, 1:(shock - 1))'; zeroResponses{shock}; ...
        info.Weights{shock}];
    [basis, upper] = qr(matrix.');
    for diagonal = (numSeries - count + 1):numSeries
        if upper(diagonal, diagonal) < 0
            basis(:, diagonal) = -basis(:, diagonal);
        end
    end
    w(offset + (1:count)) = ...
        basis(:, (numSeries - count + 1):numSeries)'*rotation(:, shock);
    offset = offset + count;
end
y = [z(1:numSeries*(info.NumCoefficients + numSeries)); w];
end

function values = localZeroRestrictions(x, info)
responses = localImpactResponses(x, info);
values = zeros(info.NumZeros, 1);
offset = 0;
for shock = 1:info.NumSeries
    count = size(responses{shock}, 1);
    values(offset + (1:count)) = responses{shock}(:, shock);
    offset = offset + count;
end
end

function responses = localImpactResponses(x, info)
numSeries = info.NumSeries;
A0 = reshape(x(1:numSeries^2), numSeries, numSeries);
impact = inv(A0).';
responses = cell(numSeries, 1);
for shock = 1:numSeries
    responses{shock} = info.Selectors{shock}*impact;
end
end

function value = localLogVolume(f, x, constraints)
derivative = localDerivative(f, x);
constraintDerivative = localDerivative(constraints, x);
nullBasis = null(constraintDerivative);
value = 0.5*localLogAbsDet((derivative*nullBasis)'*(derivative*nullBasis));
end

function derivative = localDerivative(f, x)
step = 1e-6;
baseline = f(x);
derivative = zeros(numel(baseline), numel(x));
for index = 1:numel(x)
    perturbed = x;
    perturbed(index) = perturbed(index) + step;
    derivative(:, index) = (f(perturbed) - baseline)/(2*step);
end
end

function value = localLogAbsDet(matrix)
[~, upper] = lu(matrix);
value = sum(log(abs(diag(upper))));
end
