classdef minnesotainwSpec
    %minnesotainwSpec Hyperparameter recipe for the Independent Normal-Wishart
    %   Minnesota prior (see INWBVARM; Kadiyala & Karlsson, 1997).
    %
    %   Sibling to MINNESOTASPEC, not a subclass of it - the two feed
    %   different model classes with different hyperparameter sets and,
    %   crucially, different downstream tuning stories:
    %     * MINNESOTASPEC       -> MINNESOTABVARM      -> analytic marginal
    %       likelihood -> GLPOPTIMIZEMINNESOTA (fmincon over a closed-form
    %       objective).
    %     * INWSPEC    -> INWBVARM    -> Gibbs draws, NO
    %       closed-form marginal likelihood -> hyperparameter selection needs
    %       an out-of-sample criterion (e.g. rolling-window forecast loss
    %       over a grid/search in lambda2), which is not implemented here.
    %   Sharing a base class between the two would suggest they interoperate
    %   with the same optimiser, which they do not - so they are kept
    %   separate on purpose (see the discussion that led to this file).
    %
    %   Scalar-or-bounds convention (identical to MINNESOTASPEC)
    %   -----------------------------------------------------------
    %   Each of lambda1/lambda2/lambda3 is EITHER a scalar (FIXED) or a
    %   2-element [lower upper] (FREE, to be optimised within those bounds).
    %   The validator (MUSTBESCALARORBOUNDS) is shared with MINNESOTASPEC -
    %   one implementation, not a copy.
    %
    %   No logHyperprior here (unlike MINNESOTASPEC): a hyperprior only makes
    %   sense as a correction turning a marginal-likelihood objective into a
    %   MAP objective, and this prior family has no marginal likelihood to
    %   correct.

    properties
        lambda1   (1,:) double {mustBeScalarOrBounds} = 0.2   % overall tightness
        lambda2   (1,:) double {mustBeScalarOrBounds} = 0.5   % cross-variable tightness (FREE here)
        lambda3   (1,:) double {mustBeScalarOrBounds} = 1     % lag decay
        Vc        (1,1) double {mustBePositive}      = 1e4    % constant / trend variance
        PriorMean (1,:) double = []                           % own first-lag mean ([] -> ones)
    end

    properties (Constant, Access = private)
        HyperparamNames = ["lambda1","lambda2","lambda3"]
    end

    methods

        function spec = minnesotainwSpec(nvp)
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
            %ISFREE True if the named hyperparameter is a 2-element bound.
            arguments
                spec (1,1) minnesotainwSpec
                name (1,1) string
            end
            tf = numel(spec.(name)) == 2;
        end

        function names = freeFields(spec)
            %FREEFIELDS Names of the currently-free hyperparameters, in a
            %   fixed canonical order (lambda1, lambda2, lambda3).
            arguments
                spec (1,1) minnesotainwSpec
            end
            mask  = arrayfun(@(n) spec.isFree(n), minnesotainwSpec.HyperparamNames);
            names = minnesotainwSpec.HyperparamNames(mask);
        end

        function tf = isResolved(spec)
            %ISRESOLVED True if every hyperparameter is fixed (scalar).
            arguments
                spec (1,1) minnesotainwSpec
            end
            tf = isempty(spec.freeFields());
        end

        function [x0, lb, ub, names] = pack(spec)
            %PACK Bounds and a starting point for every free hyperparameter.
            %   x0 uses the geometric mean of each [lower upper] bound - see
            %   MINNESOTASPEC.pack for the same convention.
            arguments
                spec (1,1) minnesotainwSpec
            end
            names = spec.freeFields();
            n     = numel(names);
            x0 = zeros(1, n); lb = zeros(1, n); ub = zeros(1, n);
            for i = 1:n
                b     = spec.(names(i));
                lb(i) = b(1);
                ub(i) = b(2);
                x0(i) = sqrt(b(1)*b(2));
            end
        end

        function spec = unpack(spec, x, names)
            %UNPACK Write point values back into the named fields (returns a
            %   modified COPY - value semantics, as in MINNESOTASPEC.unpack).
            arguments
                spec  (1,1) minnesotainwSpec
                x     (1,:) double
                names (1,:) string
            end
            if numel(x) ~= numel(names)
                error("minnesotainwSpec:unpack:sizeMismatch", ...
                    "x has %d elements but names has %d.", numel(x), numel(names));
            end
            for i = 1:numel(names)
                spec.(names(i)) = x(i);
            end
        end

        function mdl = build(spec, numseries, numlags, residualVariances, opts)
            %BUILD Materialise a MINNESOTAINWBVARM from this (resolved) spec.
            %   residualVariances is taken precomputed. Errors if the spec
            %   still has free fields - build needs concrete numbers, not ranges.
            arguments
                spec              (1,1) minnesotainwSpec
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
                error("minnesotainwSpec:build:notResolved", ...
                    "Cannot build: %s still free (2-element bound). Call unpack " + ...
                    "with a candidate point first.", strjoin(spec.freeFields(), ", "));
            end
            args = namedargs2cell(opts);
            mdl = minnesotainwbvarm(numseries, numlags, args{:}, ...
                ResidualVariances = residualVariances, ...
                lambda1           = spec.lambda1, ...
                lambda2           = spec.lambda2, ...
                lambda3           = spec.lambda3, ...
                Vc                = spec.Vc, ...
                PriorMean         = spec.PriorMean);
        end

    end
end
