classdef (Hidden) minnesotanbvarm < normalbvarm & minnesotabvarm & matlab.mixin.CustomDisplay
    %MINNESOTANBVARM Fixed-Sigma Normal Minnesota prior for a Bayesian VAR.
    %
    %   MINNESOTANBVARM materialises a Minnesota-shaped Normal prior with a
    %   fixed innovations covariance matrix Sigma = diag(ResidualVariances).
    %   Unlike the conjugate MNIW variant, the coefficient covariance is a full
    %   (m*n)-by-(m*n) diagonal matrix, so lambda2 can shrink cross-variable
    %   lag coefficients without affecting own lags.

    properties (SetAccess = private)
        lambda2 (1,1) double {mustBePositive} = 0.5
    end

    methods
        function obj = minnesotanbvarm(numseries, numlags, nvp, nvp2)
            arguments
                numseries (1,1) double {mustBeInteger, mustBePositive}
                numlags   (1,1) double {mustBeInteger, mustBePositive}
                nvp.ResidualVariances (1,:) double {mustBePositive} = []
                nvp.lambda1   (1,1) double {mustBePositive}    = 0.2
                nvp.lambda2   (1,1) double {mustBePositive}    = 0.5
                nvp.lambda3   (1,1) double {mustBeNonnegative} = 1
                nvp.Vc        (1,1) double {mustBePositive}    = 1e4
                nvp.PriorMean (1,:) double = []
                nvp2.Description
                nvp2.IncludeConstant
                nvp2.IncludeTrend
                nvp2.NumPredictors
                nvp2.SeriesNames
            end

            args = namedargs2cell(nvp2);
            obj@normalbvarm(numseries, numlags, args{:});

            [residualVariances, priorMean] = obj.validateMinnesotaInputs( ...
                nvp.ResidualVariances, nvp.PriorMean, "minnesotanbvarm");

            obj.ResidualVariances = residualVariances;
            obj.lambda1           = nvp.lambda1;
            obj.lambda2           = nvp.lambda2;
            obj.lambda3           = nvp.lambda3;
            obj.Vc                = nvp.Vc;
            obj.PriorMean         = priorMean;

            obj.Mu    = obj.buildMinnesotaPriorMean();
            obj.V     = obj.buildIndependentCoefficientCovariance(obj.lambda2);
            obj.Sigma = diag(obj.ResidualVariances);
        end

        function [Posterior, Summary] = estimate(obj, Y, opts)
            %ESTIMATE Analytic fixed-Sigma Normal posterior.
            arguments
                obj
                Y double {mustBeNonempty}
                opts.Display
                opts.X
                opts.Y0
            end

            args  = namedargs2cell(opts);
            [NormalPosterior, Summary] = estimate@normalbvarm(obj, Y, args{:});

            Posterior = normalbvarm(NormalPosterior.NumSeries, NormalPosterior.P, ...
                Description     = NormalPosterior.Description, ...
                SeriesNames     = NormalPosterior.SeriesNames, ...
                IncludeConstant = NormalPosterior.IncludeConstant, ...
                IncludeTrend    = NormalPosterior.IncludeTrend, ...
                NumPredictors   = NormalPosterior.NumPredictors, ...
                Mu              = NormalPosterior.Mu, ...
                V               = NormalPosterior.V, ...
                Sigma           = NormalPosterior.Sigma);
        end

        function [logML, details] = logMarginalLikelihood(obj, Y)
            %LOGMARGINALLIKELIHOOD log p(Y | hyperparameters).
            arguments
                obj
                Y double {mustBeNonempty}
            end

            [X, YResponse] = obj.regressionMatrices(Y);
            numObs = size(YResponse, 1);
            numSeries = obj.NumSeries;
            y = YResponse(:);

            design = kron(eye(numSeries), X);
            priorCovariance = (obj.V + obj.V')/2;
            innovationCovariance = kron((obj.Sigma + obj.Sigma')/2, eye(numObs));
            marginalCovariance = innovationCovariance ...
                + design*priorCovariance*design';
            marginalCovariance = (marginalCovariance + marginalCovariance')/2;

            details = struct("NumObservations", numObs);
            try
                cholFactor = chol(marginalCovariance, "lower");
                residual = y - design*obj.Mu(:);
                standardizedResidual = cholFactor \ residual;
                logDeterminant = 2*sum(log(diag(cholFactor)));

                logML = -0.5*(numel(y)*log(2*pi) ...
                    + logDeterminant ...
                    + standardizedResidual'*standardizedResidual);

                if nargout > 1
                    [posteriorMu, posteriorV] = obj.normalPosteriorMoments( ...
                        X, YResponse);
                    details = struct( ...
                        "NumObservations", numObs, ...
                        "PosteriorMu", posteriorMu, ...
                        "PosteriorV", posteriorV);
                end
            catch
                logML = -Inf;
            end
        end

    end

    methods (Access = private)
        function V = buildIndependentCoefficientCovariance(obj, lambda2)
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

        function [X, YResponse] = regressionMatrices(obj, Y)
            %REGRESSIONMATRICES VAR design matrix consistent with Mu and V.
            obj.validateData(Y);
            n      = obj.NumSeries;
            P      = obj.P;
            mm     = obj.m;
            numObs = size(Y, 1);
            T      = numObs - P;

            X = zeros(T, mm);
            for lag = 1:P
                cols = (lag - 1)*n + (1:n);
                X(:, cols) = Y((P + 1 - lag):(numObs - lag), :);
            end
            col = P*n;
            if obj.IncludeConstant
                col = col + 1;
                X(:, col) = 1;
            end
            if obj.IncludeTrend
                col = col + 1;
                X(:, col) = (1:T)';
            end
            if obj.NumPredictors > 0
                error("minnesotanbvarm:predictorsUnsupported", ...
                    ["The marginal likelihood path needs the exogenous " ...
                     "regressors, which this prior does not store. Estimate " ...
                     "with the X name-value argument instead."]);
            end

            YResponse = Y((P + 1):numObs, :);
        end

        function [posteriorMu, posteriorV] = normalPosteriorMoments(obj, X, YResponse)
            numObs = size(YResponse, 1);
            n = obj.NumSeries;
            design = kron(eye(n), X);
            y = YResponse(:);

            priorCovariance = (obj.V + obj.V')/2;
            innovationPrecision = kron(obj.Sigma \ eye(n), eye(numObs));
            priorPrecision = priorCovariance \ eye(size(priorCovariance, 1));
            posteriorPrecision = priorPrecision ...
                + design'*innovationPrecision*design;
            posteriorV = posteriorPrecision \ eye(size(posteriorPrecision, 1));
            posteriorV = (posteriorV + posteriorV')/2;
            posteriorMu = posteriorV*(priorPrecision*obj.Mu(:) ...
                + design'*innovationPrecision*y);
        end

        function validateData(obj, Y)
            if size(Y, 2) ~= obj.NumSeries
                error("minnesotanbvarm:invalidData", ...
                    "Y must have %d columns, one per series.", obj.NumSeries);
            end
            if size(Y, 1) <= obj.P
                error("minnesotanbvarm:invalidData", ...
                    "Y must have more rows than the lag order P = %d.", obj.P);
            end
        end
    end

    methods (Access = protected)
        function displayScalarObject(obj)
            disp(matlab.mixin.CustomDisplay.getSimpleHeader(obj));
            base = ["NumSeries","P","ResidualVariances","lambda1","lambda2", ...
                "lambda3","Vc","PriorMean"];
            group = obj.minnesotaPropertyGroup(base);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj, group);
        end
    end
end
