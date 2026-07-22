classdef (Abstract, Hidden) minnesotabvarmBase
    %MINNESOTABVARM Abstract shared base for Minnesota BVAR priors.
    %
    %   This class holds state and helpers common to the concrete Minnesota
    %   prior families. Construct MINNESOTAMNIWBVARM, MINNESOTAINWBVARM, or
    %   MINNESOTANBVARM instead.

    properties (SetAccess = protected)
        ResidualVariances (1,:) double {mustBePositive}
        lambda1   (1,1) double {mustBePositive}    = 0.2
        lambda3   (1,1) double {mustBeNonnegative} = 1
        Vc        (1,1) double {mustBePositive}    = 1e4
        PriorMean (1,:) double
    end

    properties (Dependent, SetAccess = private, Hidden)
        m       % total coefficients per equation = P*NumSeries + nex
        nex     % number of exogenous regressors = const + trend + predictors
    end

    methods
        function nex = get.nex(obj)
            nex = double(obj.IncludeConstant) ...
                + double(obj.IncludeTrend) ...
                + obj.NumPredictors; %#ok<*MCNPN>
        end

        function m = get.m(obj)
            m = obj.P*obj.NumSeries + obj.nex;
        end

    end

    methods (Access = protected)
        function [priorMean] = validateMinnesotaInputs( ...
                obj, residualVariancesValue, priorMeanValue, errorPrefix)
            numseries = obj.NumSeries;
 
            if numel(residualVariancesValue) ~= numseries
                error(errorPrefix + ":residualVariancesSize", ...
                    "ResidualVariances must have %d elements, one per series.", numseries);
            end

            priorMean = priorMeanValue;
            if isscalar(priorMeanValue)
                priorMean = priorMeanValue*ones(1, numseries);
            elseif numel(priorMeanValue) ~= numseries
                error(errorPrefix + ":priorMeanSize", ...
                    "PriorMean must have %d elements, one per series.", numseries);
            end

        end

        function Mu = buildMinnesotaPriorMean(obj)
            numseries     = obj.NumSeries;
            MuMat         = zeros(obj.m, numseries);
            MuMat(1:numseries, :) = diag(obj.PriorMean);
            Mu            = MuMat(:);
        end

        function group = minnesotaPropertyGroup(~, propertyNames)
            group = matlab.mixin.util.PropertyGroup(cellstr(propertyNames));
        end

        function validateMinnesotaData(obj, Y, errorPrefix)
            if size(Y, 2) ~= obj.NumSeries
                error(errorPrefix + ":invalidData", ...
                    "Y must have %d columns, one per series.", obj.NumSeries);
            end
            if size(Y, 1) <= obj.P
                error(errorPrefix + ":invalidData", ...
                    "Y must have more rows than the lag order P = %d.", obj.P);
            end
        end
    end
end
