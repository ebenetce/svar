classdef (Hidden) minnesotainwSpec < svar.minnesotaBaseSpec
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
        lambda2   (1,:) {mustBeScalarOrBounds} = 0.5   % cross-variable tightness (FREE here)
    end

    methods (Static, Hidden)
        function spec = create(varargin)
            spec = minnesotainwSpec(varargin{:});
        end
    end

    methods (Access = private)

        function spec = minnesotainwSpec(nvp)
            arguments
                nvp.lambda1   (1,:)
                nvp.lambda2   (1,:)
                nvp.lambda3   (1,:)
                nvp.Vc        (1,1) double
                nvp.PriorMean (1,:) double
            end
            spec = spec.assignSpecInputs(nvp);
        end
    end

    methods

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
            spec.assertResolvedForBuild();
            args = spec.modelConstructorArgs(residualVariances, opts);
            mdl = minnesotainwbvarm(numseries, numlags, args{:}, ...
                lambda2 = spec.lambda2);
        end

    end

    methods (Access = protected)
        function names = hyperparameterNames(~)
            names = ["lambda1","lambda2","lambda3"];
        end
    end
end
