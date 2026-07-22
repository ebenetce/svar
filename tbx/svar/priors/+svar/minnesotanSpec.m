classdef (Hidden) minnesotanSpec < svar.minnesotaBaseSpec
    %MINNESOTANSPEC Hyperparameter recipe for fixed-Sigma Normal Minnesota prior.

    properties
        lambda2   (1,:) {svar.mustBeScalarOrBounds} = 0.5
    end

    methods (Static, Hidden)
        function spec = create(varargin)
            spec = svar.minnesotanSpec(varargin{:});
        end
    end

    methods (Access = private)
        function spec = minnesotanSpec(nvp)
            arguments
                nvp.lambda1   (1,:)
                nvp.lambda2   (1,:)
                nvp.lambda3   (1,:)
                nvp.Vc        (1,1) double
                nvp.PriorMean (1,:) double = 1
            end
            spec = spec.assignSpecInputs(nvp);
        end
    end

    methods
        function mdl = build(spec, numseries, numlags, residualVariances, opts)
            arguments
                spec              (1,1) svar.minnesotanSpec
                numseries         (1,1) double {mustBeInteger, mustBePositive}
                numlags           (1,1) double {mustBeInteger, mustBePositive}
                residualVariances (1,:) double {mustBePositive}
                opts.IncludeConstant
                opts.IncludeTrend
                opts.NumPredictors
                opts.SeriesNames
                opts.Description
            end
            spec.assertResolvedForBuild();
            args = spec.modelConstructorArgs(opts);
            mdl = svar.minnesotanbvarm(numseries, numlags, residualVariances, args{:}, ...
                lambda2 = spec.lambda2);
        end
    end

    methods (Access = protected)
        function names = hyperparameterNames(~)
            names = ["lambda1","lambda2","lambda3"];
        end
    end
end
