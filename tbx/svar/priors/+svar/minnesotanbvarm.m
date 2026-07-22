classdef minnesotanbvarm < normalbvarm & svar.minnesotabvarmBase & matlab.mixin.CustomDisplay
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
        function obj = minnesotanbvarm(numseries, numlags, residualVariances, nvp, nvp2)
            arguments
                numseries (1,1) double {mustBeInteger, mustBePositive}
                numlags   (1,1) double {mustBeInteger, mustBePositive}
                residualVariances (1,:) double {mustBePositive}
                nvp.lambda1   (1,1) double {mustBePositive}    = 0.2
                nvp.lambda2   (1,1) double {mustBePositive}    = 0.5
                nvp.lambda3   (1,1) double {mustBeNonnegative} = 1
                nvp.Vc        (1,1) double {mustBePositive}    = 1e4
                nvp.PriorMean (1,:) double = ones(1, numseries);
                nvp2.Description
                nvp2.IncludeConstant
                nvp2.IncludeTrend
                nvp2.NumPredictors
                nvp2.SeriesNames
            end

            args = namedargs2cell(nvp2);
            obj@normalbvarm(numseries, numlags, args{:});

            priorMean = obj.validateMinnesotaInputs( ...
                residualVariances, nvp.PriorMean, "minnesotanbvarm");

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
