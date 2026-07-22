classdef minnesotamniwSpec < svar.minnesotaBaseSpec
    %MINNESOTASPEC Hyperparameter recipe for a conjugate Minnesota BVAR prior.
    %
    %   A lightweight value object that holds the Minnesota hyperparameters
    %   and knows how to (a) materialise a concrete MINNESOTABVARM from itself
    %   plus a residual-variance vector, and (b) pack/unpack its FREE fields
    %   to/from a plain bounded vector for FMINCON.
    %
    %   Scalar-or-bounds convention
    %   ----------------------------
    %   Each of lambda1/lambda3/lambda4/lambda5 is EITHER:
    %     * a scalar  -> FIXED at that value (not optimised), or
    %     * a 2-element [lower upper] -> FREE, optimised within those bounds.
    %   This mirrors the classic GLP-style hyperparameter-search interface:
    %   pass a number to hold it fixed, pass a range to tune it. Which fields
    %   are free is therefore implicit in the spec itself - no separate
    %   FreeParams list to keep in sync.
    %
    %   Separation of concerns
    %   ----------------------
    %   * The spec is the mutable RECIPE: hyperparameter vocabulary, defaults,
    %     the free/fixed convention above, and the build step. It holds no
    %     data and no moments, and it never carries a per-series psi (psi's
    %     dimension depends on NumSeries, which the spec deliberately does not
    %     know - see the GLP wrapper for how a psi band is optimised
    %     alongside the spec's free lambdas).
    %   * MINNESOTABVARM is the materialised RESULT. build() copies the
    %     hyperparameter VALUES into it, so mutating the spec afterwards never
    %     disturbs an already-built model - no shared state.
    %   * logHyperprior is a method here (not on minnesotabvarm): it is a
    %     belief about the lambdas themselves, needs only the spec's scalar
    %     values, and requires no build step to evaluate.
    %   * A GLP-style optimiser is the third layer (see GLPOPTIMIZEMINNESOTA):
    %     it reads the spec's free fields and bounds via pack(), lets FMINCON
    %     search them (plus, separately, a psi band), and rebuilds via
    %     build()/unpack() at each candidate point.
    %
    %   Off switches: lambda4 = lambda5 = Inf (as a FIXED scalar) disables the
    %   sum-of-coefficients and dummy-initial-observation priors respectively
    %   (see MINNESOTABVARM). Inf is not a valid bound endpoint.
    %
    %   Note on lambda2: there is no separate cross-variable tightness field.
    %   Under the shared-diagonal-V / diag(psi) conjugate structure this class
    %   uses - algebraically identical to the dummy-observation construction
    %   this design was checked against - the implied coefficient variance is
    %   lambda1^2/l^(2*lambda3) (own) and lambda1^2/l^(2*lambda3) * (psi_j/psi_i)
    %   (cross): both governed by the SAME lambda1/lambda3, with only the
    %   psi ratio distinguishing them. There is no free parameter left over to
    %   call lambda2; it is not omitted for convenience, it does not exist
    %   under this structure.

    properties
        lambda4   (1,:) {svar.mustBeScalarOrBounds} = Inf   % sum-of-coeff (Inf = off)
        lambda5   (1,:) {svar.mustBeScalarOrBounds} = Inf   % dummy-init-obs (Inf = off)
    end

    methods (Static, Hidden)
        function spec = create(varargin)
            spec = svar.minnesotamniwSpec(varargin{:});
        end
    end

    methods (Access = private)

        function spec = minnesotamniwSpec(nvp)
            arguments
                nvp.lambda1   (1,:)
                nvp.lambda3   (1,:)
                nvp.lambda4   (1,:)
                nvp.lambda5   (1,:)
                nvp.Vc        (1,1) double
                nvp.PriorMean (1,:) double = 1
            end
            spec = spec.assignSpecInputs(nvp);
        end
    end

    methods
        function mdl = build(spec, numseries, numlags, residualVariances, opts)
            %BUILD Materialise a MINNESOTABVARM from this (resolved) spec.
            %   residualVariances is taken precomputed (compute it once,
            %   outside any tuning loop). Errors if the spec still has free
            %   fields - build needs concrete numbers, not ranges.
            arguments
                spec              (1,1) svar.minnesotamniwSpec
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
            mdl = svar.minnesotamniwbvarm(numseries, numlags, residualVariances, args{:}, ...
                lambda4 = spec.lambda4, ...
                lambda5 = spec.lambda5);
        end

        function lp = logHyperprior(spec, x, names)
            %LOGHYPERPRIOR Log density from embedded lambda hyperpriors.
            arguments
                spec  (1,1) svar.minnesotamniwSpec
                x     (1,:) double {mustBeFinite}
                names (1,:) string
            end

            lp = logHyperprior@svar.minnesotaBaseSpec(spec, x, names);
        end

    end

    methods (Access = protected)
        function names = hyperparameterNames(~)
            names = ["lambda1","lambda3","lambda4","lambda5"];
        end
    end
end
