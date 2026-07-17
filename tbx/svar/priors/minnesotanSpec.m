classdef minnesotanSpec
    %MINNESOTANSPEC Hyperparameter recipe for fixed-Sigma Normal Minnesota prior.

    properties
        lambda1   (1,:) double {mustBeScalarOrBounds} = 0.2
        lambda2   (1,:) double {mustBeScalarOrBounds} = 0.5
        lambda3   (1,:) double {mustBeScalarOrBounds} = 1
        Vc        (1,1) double {mustBePositive}       = 1e4
        PriorMean (1,:) double = []
    end

    properties (Constant, Access = private)
        HyperparamNames = ["lambda1","lambda2","lambda3"]
    end

    methods
        function spec = minnesotanSpec(nvp)
            arguments
                nvp.lambda1   (1,:) double
                nvp.lambda2   (1,:) double
                nvp.lambda3   (1,:) double
                nvp.Vc        (1,1) double
                nvp.PriorMean (1,:) double
            end
            for f = string(fieldnames(nvp))'
                spec.(f) = nvp.(f);
            end
        end

        function tf = isFree(spec, name)
            arguments
                spec (1,1) minnesotanSpec
                name (1,1) string
            end
            tf = numel(spec.(name)) == 2;
        end

        function names = freeFields(spec)
            arguments
                spec (1,1) minnesotanSpec
            end
            mask  = arrayfun(@(n) spec.isFree(n), minnesotanSpec.HyperparamNames);
            names = minnesotanSpec.HyperparamNames(mask);
        end

        function tf = isResolved(spec)
            arguments
                spec (1,1) minnesotanSpec
            end
            tf = isempty(spec.freeFields());
        end

        function [x0, lb, ub, names] = pack(spec)
            arguments
                spec (1,1) minnesotanSpec
            end
            names = spec.freeFields();
            n     = numel(names);
            x0 = zeros(1, n);
            lb = zeros(1, n);
            ub = zeros(1, n);
            for i = 1:n
                b     = spec.(names(i));
                lb(i) = b(1);
                ub(i) = b(2);
                x0(i) = sqrt(b(1)*b(2));
            end
        end

        function spec = unpack(spec, x, names)
            arguments
                spec  (1,1) minnesotanSpec
                x     (1,:) double
                names (1,:) string
            end
            if numel(x) ~= numel(names)
                error("minnesotanSpec:unpack:sizeMismatch", ...
                    "x has %d elements but names has %d.", numel(x), numel(names));
            end
            for i = 1:numel(names)
                spec.(names(i)) = x(i);
            end
        end

        function mdl = build(spec, numseries, numlags, residualVariances, opts)
            arguments
                spec              (1,1) minnesotanSpec
                numseries         (1,1) double {mustBeInteger, mustBePositive}
                numlags           (1,1) double {mustBeInteger, mustBePositive}
                residualVariances (1,:) double {mustBePositive}
                opts.IncludeConstant
                opts.IncludeTrend
                opts.NumPredictors
                opts.SeriesNames
                opts.Description
            end
            if ~spec.isResolved()
                error("minnesotanSpec:build:notResolved", ...
                    "Cannot build: %s still free (2-element bound). Call unpack " + ...
                    "with a candidate point first.", strjoin(spec.freeFields(), ", "));
            end
            args = namedargs2cell(opts);
            mdl = minnesotanbvarm(numseries, numlags, args{:}, ...
                ResidualVariances = residualVariances, ...
                lambda1           = spec.lambda1, ...
                lambda2           = spec.lambda2, ...
                lambda3           = spec.lambda3, ...
                Vc                = spec.Vc, ...
                PriorMean         = spec.PriorMean);
        end
    end
end
