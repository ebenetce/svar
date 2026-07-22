classdef (Abstract) minnesotaBaseSpec
    %MINNESOTABASESPEC Shared scalar-or-bounds behavior for Minnesota specs.

    properties
        lambda1   (1,:) {svar.mustBeScalarOrBounds} = 0.2
        lambda3   (1,:) {svar.mustBeScalarOrBounds} = 1
        Vc        (1,1) double {mustBePositive}       = 1e4
        PriorMean (1,:) double = []
    end

    methods
        function tf = isFree(spec, name)
            arguments
                spec (1,1) svar.minnesotaBaseSpec
                name (1,1) string
            end
            tf = numel(spec.(name)) == 2 || isa(spec.(name), "hyperprior");
        end

        function names = freeFields(spec)
            arguments
                spec (1,1) svar.minnesotaBaseSpec
            end
            hyperparameters = spec.hyperparameterNames();
            mask = arrayfun(@(n) spec.isFree(n), hyperparameters);
            names = hyperparameters(mask);
        end

        function tf = isResolved(spec)
            arguments
                spec (1,1) svar.minnesotaBaseSpec
            end
            tf = isempty(spec.freeFields());
        end

        function [x0, lb, ub, names] = pack(spec)
            arguments
                spec (1,1) svar.minnesotaBaseSpec
            end
            names = spec.freeFields();
            n = numel(names);
            x0 = zeros(1, n);
            lb = zeros(1, n);
            ub = zeros(1, n);
            for i = 1:n
                b = spec.(names(i));
                if isa(b, 'hyperprior')
                    lb(i) = b.Bounds(1);
                    ub(i) = b.Bounds(2);
                    x0(i) = b.X0;
                else
                    lb(i) = b(1);
                    ub(i) = b(2);
                    x0(i) = sqrt(b(1)*b(2));
                end
            end
        end

        function spec = unpack(spec, x, names)
            arguments
                spec  (1,1) svar.minnesotaBaseSpec
                x     (1,:) double
                names (1,:) string
            end
            if numel(x) ~= numel(names)
                error(string(class(spec)) + ":unpack:sizeMismatch", ...
                    "x has %d elements but names has %d.", numel(x), numel(names));
            end
            for i = 1:numel(names)
                spec.(names(i)) = x(i);
            end
        end

        function lp = logHyperprior(spec, x, names)
            arguments
                spec  (1,1) svar.minnesotaBaseSpec
                x     (1,:) double {mustBeFinite}
                names (1,:) string
            end
            if numel(x) ~= numel(names)
                error(string(class(spec)) + ":logHyperprior:sizeMismatch", ...
                    "x has %d elements but names has %d.", numel(x), numel(names));
            end
            lp = 0;
            for i = 1:numel(names)
                prior = spec.(names(i));
                if isa(prior, "hyperprior")
                    lp = lp + prior.logpdf(x(i));
                end
            end
        end

        function names = hyperpriorFields(spec)
            arguments
                spec (1,1) svar.minnesotaBaseSpec
            end
            hyperparameters = spec.hyperparameterNames();
            mask = arrayfun(@(n) isa(spec.(n), "hyperprior"), hyperparameters);
            names = hyperparameters(mask);
        end
    end

    methods (Access = protected)
        function spec = assignSpecInputs(spec, inputs)
            for f = string(fieldnames(inputs))'
                spec.(f) = inputs.(f);
            end
        end

        function assertResolvedForBuild(spec)
            if ~spec.isResolved()
                error(string(class(spec)) + ":build:notResolved", ...
                    "Cannot build: %s still free. Call unpack " + ...
                    "with a candidate point first.", strjoin(spec.freeFields(), ", "));
            end
        end

        function args = modelConstructorArgs(spec, opts)
            arguments
                spec (1,1) svar.minnesotaBaseSpec
                opts (1,1) struct
            end
            args = [namedargs2cell(opts), ...
                {"lambda1", spec.lambda1, ...
                "lambda3", spec.lambda3, ...
                "Vc", spec.Vc, ...
                "PriorMean", spec.PriorMean}];
        end
    end

    methods (Abstract, Access = protected)
        names = hyperparameterNames(spec)
    end
end
