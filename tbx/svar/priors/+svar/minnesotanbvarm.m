classdef minnesotanbvarm < normalbvarm & svar.minnesotabvarmBase & matlab.mixin.CustomDisplay
    %MINNESOTANBVARM Fixed-Sigma Normal Minnesota prior for a Bayesian VAR.
    %
    %   PriorMdl = SVAR.MINNESOTANBVARM(NUMSERIES,NUMLAGS,Y) materialises a
    %   Minnesota-shaped Normal prior for the sample Y, with a FIXED
    %   innovations covariance Sigma = diag(Psi). Prefer the MINNESOTABVARM
    %   front door, which reaches this class as Type="normal".
    %
    %   Like the INW family and unlike the conjugate one, the coefficient
    %   covariance is a full (m*n)-by-(m*n) diagonal matrix, so lambda2 can
    %   shrink cross-variable lag coefficients without affecting own lags.
    %   Sigma is not estimated - it is pinned at the residual-variance scale -
    %   which is what makes the posterior analytic despite the non-Kronecker V.
    %
    %   See also MINNESOTABVARM, SVAR.MINNESOTAINWBVARM, NORMALBVARM.

    properties (SetAccess = private)
        lambda2 (1,1) double {mustBePositive} = 0.5
    end

    methods
        function obj = minnesotanbvarm(numseries, numlags, Y, nvp, nvp2)
            arguments
                numseries (1,1) double {mustBeInteger, mustBePositive}
                numlags   (1,1) double {mustBeInteger, mustBePositive}
                Y                {mustBeNonempty}
                nvp.Psi           = "exact"
                nvp.lambda1   (1,1) double {mustBePositive}    = 0.2
                nvp.lambda2   (1,1) double {mustBePositive}    = 0.5
                nvp.lambda3   (1,1) double {mustBeNonnegative} = 1
                nvp.Vc        (1,1) double {mustBePositive}    = 1e4
                nvp.PriorMean (1,:) double = ones(1, numseries)
                nvp2.Description
                nvp2.IncludeConstant
                nvp2.IncludeTrend
                nvp2.NumPredictors
                nvp2.SeriesNames
            end

            [Y, nvp2] = svar.minnesotanbvarm.resolveSample( ...
                Y, numseries, numlags, nvp2, "minnesotanbvarm");

            args = namedargs2cell(nvp2);
            obj@normalbvarm(numseries, numlags, args{:});

            obj = obj.configureMinnesota(Y, numlags, nvp, "minnesotanbvarm");
            obj.lambda2 = nvp.lambda2;

            obj.Mu    = obj.buildMinnesotaPriorMean();
            obj.V     = obj.buildIndependentCoefficientCovariance(obj.lambda2);
            obj.Sigma = diag(obj.ResidualVariances);
        end

        function [Posterior, Summary] = estimate(obj, varargin)
            %ESTIMATE Analytic fixed-Sigma Normal posterior for the stored sample.
            [Y, args] = obj.peelSample(varargin, "estimate", "minnesotanbvarm");
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

        function varargout = simulate(obj, varargin)
            %SIMULATE Draw coefficients given the stored sample.
            [Y, args] = obj.peelSample(varargin, "simulate", "minnesotanbvarm");
            [varargout{1:nargout}] = simulate@normalbvarm(obj, Y, args{:});
        end

        function varargout = forecast(obj, numperiods, varargin)
            %FORECAST Forecast responses beyond the stored sample.
            arguments
                obj        (1,1) svar.minnesotanbvarm
                numperiods (1,1) double {mustBeInteger, mustBePositive}
            end
            arguments (Repeating)
                varargin
            end
            [Y, args] = obj.peelSample(varargin, "forecast", "minnesotanbvarm");
            [varargout{1:nargout}] = ...
                forecast@normalbvarm(obj, numperiods, Y, args{:});
        end

        function varargout = simsmooth(obj, varargin)
            %SIMSMOOTH Simulation smoother, defaulting to the stored sample.
            %   Not sample-checked - see SVAR.MINNESOTAMNIWBVARM/SIMSMOOTH.
            [Y, args] = obj.peelSample(varargin, "simsmooth", ...
                "minnesotanbvarm", false);
            [varargout{1:nargout}] = simsmooth@normalbvarm(obj, Y, args{:});
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
