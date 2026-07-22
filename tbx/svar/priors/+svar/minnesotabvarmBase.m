classdef (Abstract) minnesotabvarmBase
    %MINNESOTABVARMBASE Abstract shared base for Minnesota BVAR priors.
    %
    %   Holds the state and helpers common to the concrete Minnesota prior
    %   families. Construct through MINNESOTABVARM, or call
    %   SVAR.MINNESOTAMNIWBVARM / SVAR.MINNESOTAINWBVARM / SVAR.MINNESOTANBVARM
    %   directly.
    %
    %   The estimation sample lives here
    %   -------------------------------
    %   Every Minnesota family needs data before it is a usable prior: the
    %   residual-variance scale Psi is estimated from Y in all three, and the
    %   conjugate family additionally derives its lambda4/lambda5 dummy
    %   observations from it. So Y is taken by the constructor and retained
    %   here, once, rather than by each subclass - which also puts the
    %   sample-resolution, Psi-resolution, and sample-consistency helpers in
    %   one place instead of three.
    %
    %   Because the sample is stored, ESTIMATE/SIMULATE/FORECAST and
    %   LOGMARGINALLIKELIHOOD may be called with no data argument. Passing a
    %   DIFFERENT sample is an error: Psi and any dummy observations derive
    %   from the stored one, so a different sample needs a different object.
    %
    %   See also MINNESOTABVARM, ESTIMATERESIDUALVARIANCES.

    properties (SetAccess = protected)
        Y                 double                     % estimation sample
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

    methods (Static, Access = protected)

        function [Y, nvp2] = resolveSample(Y, numseries, numlags, nvp2, errorPrefix)
            %RESOLVESAMPLE Validate Y and adopt tabular variable names.
            %   Called BEFORE the superclass constructor, hence static, taking
            %   the shape parameters explicitly rather than reading them off
            %   an object that does not exist yet.
            if istabular(Y)
                if ~isfield(nvp2, "SeriesNames")
                    nvp2.SeriesNames = string(Y.Properties.VariableNames);
                end
                Y = Y{:,:};
            end

            if ~isnumeric(Y) || ~ismatrix(Y)
                error(errorPrefix + ":invalidData", ...
                    "Y must be a numeric matrix or a table/timetable.");
            end
            if size(Y, 2) ~= numseries
                error(errorPrefix + ":invalidData", ...
                    "Y must have %d columns, one per series.", numseries);
            end
            if size(Y, 1) <= numlags
                error(errorPrefix + ":invalidData", ...
                    "Y must have more rows than the lag order P = %d.", numlags);
            end
            Y = double(Y);
        end

    end

    methods (Access = protected)

        function obj = configureMinnesota(obj, Y, numlags, nvp, errorPrefix)
            %CONFIGUREMINNESOTA Resolve and store the shared Minnesota state.
            %   Subclasses call this straight after their superclass
            %   constructor, then assign their own hyperparameters and build
            %   their own moments.
            residualVariances = obj.resolvePsi(nvp.Psi, Y, numlags, errorPrefix);
            priorMean = obj.validateMinnesotaInputs( ...
                residualVariances, nvp.PriorMean, errorPrefix);

            obj.Y                 = Y;
            obj.ResidualVariances = residualVariances;
            obj.lambda1           = nvp.lambda1;
            obj.lambda3           = nvp.lambda3;
            obj.Vc                = nvp.Vc;
            obj.PriorMean         = priorMean;
        end

        function psi = resolvePsi(obj, Psi, Y, numlags, errorPrefix)
            %RESOLVEPSI Residual-variance scale from a vector or an estimator name.
            if isstring(Psi) || ischar(Psi)
                method = string(Psi);
                mustBeMember(method, ["exact","conditional"]);
                psi = estimateResidualVariances(Y, numlags, Method=method);
                return
            end

            if ~isnumeric(Psi) || ~isvector(Psi) || numel(Psi) ~= obj.NumSeries
                error(errorPrefix + ":invalidPsi", ...
                    "Psi must be a 1-by-%d numeric vector or one of " + ...
                    """exact"" / ""conditional"".", obj.NumSeries);
            end
            if any(~isfinite(Psi)) || any(Psi <= 0)
                error(errorPrefix + ":invalidPsi", ...
                    "Psi values must be finite and positive.");
            end
            psi = reshape(double(Psi), 1, []);
        end

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

        function V = buildIndependentCoefficientCovariance(obj, lambda2)
            %BUILDINDEPENDENTCOEFFICIENTCOVARIANCE Full (m*n)-by-(m*n) diagonal V.
            %   Shared by the INW and fixed-Sigma Normal families: both drop
            %   the Kronecker link to Omega, so V is indexed directly over
            %   (regressor, target) PAIRS and lambda2 can shrink cross-variable
            %   lags without touching own lags.
            %
            %       Var(B_{lag l, source k -> target i}) =
            %           lambda1^2 / l^(2*lambda3) * (psi_i/psi_k)              own
            %           lambda1^2 * lambda2^2 / l^(2*lambda3) * (psi_i/psi_k)  cross
            n   = obj.NumSeries;
            P   = obj.P;
            mm  = obj.m;
            psi = obj.ResidualVariances(:);

            vDiag = zeros(mm*n, 1);
            for i = 1:n
                for r = 1:mm
                    if r <= P*n
                        lag    = ceil(r / n);
                        source = mod(r - 1, n) + 1;
                        crossFactor = 1;
                        if source ~= i
                            crossFactor = lambda2^2;
                        end
                        variance = crossFactor * obj.lambda1^2 ...
                            / lag^(2*obj.lambda3) * (psi(i) / psi(source));
                    else
                        variance = obj.Vc;
                    end
                    vDiag((i - 1)*mm + r) = variance;
                end
            end

            V = diag(vDiag);
        end

        function group = minnesotaPropertyGroup(~, propertyNames)
            group = matlab.mixin.util.PropertyGroup(cellstr(propertyNames));
        end

    end
end
