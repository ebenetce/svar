classdef minnesotamniwSpec
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
        lambda1   (1,:) double {mustBeScalarOrBounds} = 0.2   % overall tightness
        lambda3   (1,:) double {mustBeScalarOrBounds} = 1     % lag decay
        lambda4   (1,:) double {mustBeScalarOrBounds} = Inf   % sum-of-coeff (Inf = off)
        lambda5   (1,:) double {mustBeScalarOrBounds} = Inf   % dummy-init-obs (Inf = off)
        Vc        (1,1) double {mustBePositive}      = 1e4    % constant / trend variance
        PriorMean (1,:) double = []                           % own first-lag mean ([] -> ones)
    end

    properties (Constant, Access = private)
        HyperparamNames = ["lambda1","lambda3","lambda4","lambda5"]
    end

    methods

        function spec = minnesotamniwSpec(nvp)
            arguments
                nvp.lambda1   (1,:) double
                nvp.lambda3   (1,:) double
                nvp.lambda4   (1,:) double
                nvp.lambda5   (1,:) double
                nvp.Vc        (1,1) double
                nvp.PriorMean (1,:) double
            end
            % Assign only the fields the caller actually provided; the rest
            % keep their property defaults above. Property validators enforce
            % the scalar-or-bounds convention on assignment.
            for f = string(fieldnames(nvp))'
                spec.(f) = nvp.(f);
            end
        end

        function tf = isFree(spec, name)
            %ISFREE True if the named hyperparameter is a 2-element bound.
            arguments
                spec (1,1) minnesotaSpec
                name (1,1) string
            end
            tf = numel(spec.(name)) == 2;
        end

        function names = freeFields(spec)
            %FREEFIELDS Names of the currently-free hyperparameters, in a
            %   fixed canonical order (lambda1, lambda3, lambda4, lambda5).
            arguments
                spec (1,1) minnesotaSpec
            end
            mask  = arrayfun(@(n) spec.isFree(n), minnesotamniwSpec.HyperparamNames);
            names = minnesotamniwSpec.HyperparamNames(mask);
        end

        function tf = isResolved(spec)
            %ISRESOLVED True if every hyperparameter is fixed (scalar).
            %   build() and logHyperprior() require a resolved spec - neither
            %   can evaluate against a range.
            arguments
                spec (1,1) minnesotaSpec
            end
            tf = isempty(spec.freeFields());
        end

        function [x0, lb, ub, names] = pack(spec)
            %PACK Bounds and a starting point for every free hyperparameter.
            %   x0 uses the geometric mean of each [lower upper] bound, the
            %   natural default for a positive, typically-log-scaled
            %   hyperparameter (e.g. a symmetric multiplicative band around 1
            %   has geometric-mean midpoint exactly 1).
            arguments
                spec (1,1) minnesotaSpec
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
            %UNPACK Write point values back into the named fields.
            %   Returns a modified COPY (value semantics), leaving the
            %   caller's spec untouched. Each field becomes FIXED (scalar) at
            %   the supplied value, regardless of whether it was free before -
            %   this is how a candidate point during optimisation becomes a
            %   concrete, buildable spec.
            arguments
                spec  (1,1) minnesotaSpec
                x     (1,:) double
                names (1,:) string
            end
            if numel(x) ~= numel(names)
                error("minnesotaSpec:unpack:sizeMismatch", ...
                    "x has %d elements but names has %d.", numel(x), numel(names));
            end
            for i = 1:numel(names)
                spec.(names(i)) = x(i);
            end
        end

        function mdl = build(spec, numseries, numlags, ppsi, opts)
            %BUILD Materialise a MINNESOTABVARM from this (resolved) spec.
            %   ppsi is taken precomputed (compute it once, outside any
            %   tuning loop). Errors if the spec still has free fields - build
            %   needs concrete numbers, not ranges.
            arguments
                spec       (1,1) minnesotaSpec
                numseries  (1,1) double {mustBeInteger, mustBePositive}
                numlags    (1,1) double {mustBeInteger, mustBePositive}
                ppsi       (1,:) double {mustBePositive}
                opts.IncludeConstant
                opts.IncludeTrend
                opts.NumPredictors
                opts.SeriesNames
                opts.Description
            end
            if ~spec.isResolved()
                error("minnesotaSpec:build:notResolved", ...
                    "Cannot build: %s still free (2-element bound). Call unpack " + ...
                    "with a candidate point first.", strjoin(spec.freeFields(), ", "));
            end
            args = namedargs2cell(opts);
            mdl = minnesotamniwbvarm(numseries, numlags, args{:}, ...
                ppsi      = ppsi, ...
                lambda1   = spec.lambda1, ...
                lambda3   = spec.lambda3, ...
                lambda4   = spec.lambda4, ...
                lambda5   = spec.lambda5, ...
                Vc        = spec.Vc, ...
                PriorMean = spec.PriorMean);
        end

        function lp = logHyperprior(spec, priorcoef, ppsi)
            %LOGHYPERPRIOR Log prior density on the Minnesota hyperparameters.
            %   lp = spec.logHyperprior(priorcoef)        % lambdas only
            %   lp = spec.logHyperprior(priorcoef, ppsi)  % lambdas + psi
            %
            %   Returns the POSITIVE log density log p(lambda1, lambda4,
            %   lambda5, [psi]). Needs only the spec's scalar hyperparameter
            %   values - no build step, no data - which is why it lives here
            %   rather than on minnesotabvarm. Errors if the spec still has
            %   free fields (a density needs a point, not a range).
            %
            %   This is the piece that turns a marginal-likelihood objective
            %   into a MAP (posterior-mode) objective. Compose them in the
            %   optimiser:
            %
            %       mdl    = spec.build(n, p, psi, ...);
            %       negObj = mdl.negativeLogMarginalLikelihood(Y) ...
            %                - spec.logHyperprior(priorcoef, psi);
            %
            %   Omit the logHyperprior term for pure marginal-likelihood tuning.
            %
            %   priorcoef is a struct whose OPTIONAL fields switch each term on:
            %       priorcoef.lambda1.k , priorcoef.lambda1.theta   % Gamma(shape,scale)
            %       priorcoef.lambda4.k , priorcoef.lambda4.theta   % Gamma(shape,scale)
            %       priorcoef.lambda5.k , priorcoef.lambda5.theta   % Gamma(shape,scale)
            %       priorcoef.psi.alpha , priorcoef.psi.beta        % Inverse-Gamma(shape,scale)
            %   An absent field contributes nothing (flat prior on that
            %   parameter). A dummy switched off (Inf) contributes nothing even
            %   if its priorcoef field is present.
            arguments
                spec      (1,1) minnesotaSpec
                priorcoef (1,1) struct
                ppsi      (1,:) double = []
            end

            if ~spec.isResolved()
                error("minnesotaSpec:logHyperprior:notResolved", ...
                    "Cannot evaluate a density on a range: %s still free.", ...
                    strjoin(spec.freeFields(), ", "));
            end

            lp = 0;

            if isfield(priorcoef, "lambda1")
                lp = lp + minnesotamniwSpec.gammaLogPdf( ...
                    spec.lambda1, priorcoef.lambda1.k, priorcoef.lambda1.theta);
            end

            if isfield(priorcoef, "lambda4") && isfinite(spec.lambda4)
                lp = lp + minnesotamniwSpec.gammaLogPdf( ...
                    spec.lambda4, priorcoef.lambda4.k, priorcoef.lambda4.theta);
            end

            if isfield(priorcoef, "lambda5") && isfinite(spec.lambda5)
                lp = lp + minnesotamniwSpec.gammaLogPdf( ...
                    spec.lambda5, priorcoef.lambda5.k, priorcoef.lambda5.theta);
            end

            if isfield(priorcoef, "psi") && ~isempty(ppsi)
                lp = lp + sum(minnesotamniwSpec.invGammaLogPdf( ...
                    ppsi, priorcoef.psi.alpha, priorcoef.psi.beta));
            end
        end

    end

    % ---- density helpers (no state, kept private to this class) ---------
    methods (Static, Access = private)

        function r = gammaLogPdf(x, k, theta)
            %GAMMALOGPDF Log density of Gamma(shape k, scale theta) at x > 0.
            r = (k - 1).*log(x) - x./theta - k.*log(theta) - gammaln(k);
        end

        function r = invGammaLogPdf(x, alpha, beta)
            %INVGAMMALOGPDF Log density of Inverse-Gamma(shape alpha, scale beta) at x > 0.
            r = alpha.*log(beta) - (alpha + 1).*log(x) - beta./x - gammaln(alpha);
        end

    end
end

% mustBeScalarOrBounds is defined in its own file (mustBeScalarOrBounds.m),
% shared with MINNESOTAINWSPEC - see that file for the validator itself.