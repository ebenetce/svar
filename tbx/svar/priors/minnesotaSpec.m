classdef minnesotaSpec
    %MINNESOTASPEC Hyperparameter recipe for a conjugate Minnesota BVAR prior.
    %
    %   A lightweight value object that holds the Minnesota hyperparameters
    %   and knows how to (a) materialise a concrete MINNESOTABVARM from itself
    %   plus a residual-variance vector, and (b) pack/unpack a subset of its
    %   fields to/from an unconstrained optimiser vector.
    %
    %   Separation of concerns
    %   ----------------------
    %   * The spec is the mutable RECIPE. It owns the hyperparameter vocabulary,
    %     the defaults, the positivity reparametrisation (log-space, in pack/
    %     unpack), and the build step. It holds no data and no moments.
    %   * MINNESOTABVARM is the materialised RESULT. build() copies the
    %     hyperparameter VALUES into it (value semantics), so mutating the spec
    %     afterwards never disturbs an already-built model - no shared state.
    %   * logHyperprior is also a method here, not on minnesotabvarm: it is a
    %     belief about the lambdas themselves (independent of any data or
    %     built model), needs only the spec's scalar values, and requires no
    %     build step to evaluate - so it belongs with pack/unpack, not bolted
    %     onto the materialised prior.
    %   * A GLP-style optimiser is the third layer: it varies a subset of the
    %     spec's fields (pack/unpack), rebuilds cheaply for each candidate, and
    %     evaluates minnesotabvarm.negativeLogMarginalLikelihood - optionally
    %     minus spec.logHyperprior for a MAP objective. The optimiser owns none
    %     of the algebra; it only chooses lambdas.
    %
    %   Off switches: lambda4 = lambda5 = Inf disables the sum-of-coefficients
    %   and dummy-initial-observation priors respectively (see MINNESOTABVARM).

    properties
        lambda1   (1,1) double {mustBePositive}    = 0.2   % overall tightness
        lambda3   (1,1) double {mustBeNonnegative} = 1     % lag decay
        lambda4   (1,1) double {mustBePositive}    = Inf   % sum-of-coeff (Inf = off)
        lambda5   (1,1) double {mustBePositive}    = Inf   % dummy-init-obs (Inf = off)
        Vc        (1,1) double {mustBePositive}    = 1e4   % constant / trend variance
        PriorMean (1,:) double = []                        % own first-lag mean ([] -> ones)
    end

    methods

        function spec = minnesotaSpec(nvp)
            arguments
                nvp.lambda1   (1,1) double {mustBePositive}
                nvp.lambda3   (1,1) double {mustBeNonnegative}
                nvp.lambda4   (1,1) double {mustBePositive}
                nvp.lambda5   (1,1) double {mustBePositive}
                nvp.Vc        (1,1) double {mustBePositive}
                nvp.PriorMean (1,:) double
            end
            % Assign only the fields the caller actually provided; the rest
            % keep their property defaults above.
            for f = string(fieldnames(nvp))'
                spec.(f) = nvp.(f);
            end
        end

        function mdl = build(spec, numseries, numlags, ppsi, opts)
            %BUILD Materialise a MINNESOTABVARM from this spec and a ppsi vector.
            %   ppsi is taken precomputed (compute it once, outside any tuning
            %   loop). opts pass through to the conjugate superclass layout.
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
            args = namedargs2cell(opts);
            mdl = minnesotabvarm(numseries, numlags, ...
                args{:}, ...
                ppsi      = ppsi, ...
                lambda1   = spec.lambda1, ...
                lambda3   = spec.lambda3, ...
                lambda4   = spec.lambda4, ...
                lambda5   = spec.lambda5, ...
                Vc        = spec.Vc, ...
                PriorMean = spec.PriorMean);
        end

        function theta = pack(spec, freeParams)
            %PACK Map the named free hyperparameters to an unconstrained vector.
            %   Log-space keeps them positive for an unconstrained optimiser
            %   (fminsearch/fminunc). Only pass ACTIVE (finite) parameters as
            %   free - an Inf here would produce Inf in theta.
            arguments
                spec       (1,1) minnesotaSpec
                freeParams (1,:) string
            end
            theta = zeros(1, numel(freeParams));
            for i = 1:numel(freeParams)
                value = spec.(freeParams(i));
                if ~isfinite(value) || value <= 0
                    error("minnesotaSpec:pack:notOptimisable", ...
                        "Free parameter '%s' must be finite and positive; got %g.", ...
                        freeParams(i), value);
                end
                theta(i) = log(value);
            end
        end

        function spec = unpack(spec, theta, freeParams)
            %UNPACK Write an optimiser vector back into the named free fields.
            %   Returns a modified COPY (value semantics), leaving the caller's
            %   spec untouched - so each optimiser probe is independent.
            arguments
                spec       (1,1) minnesotaSpec
                theta      (1,:) double
                freeParams (1,:) string
            end
            if numel(theta) ~= numel(freeParams)
                error("minnesotaSpec:unpack:sizeMismatch", ...
                    "theta has %d elements but freeParams has %d.", ...
                    numel(theta), numel(freeParams));
            end
            for i = 1:numel(freeParams)
                spec.(freeParams(i)) = exp(theta(i));
            end
        end

        function lp = logHyperprior(spec, priorcoef, ppsi)
            %LOGHYPERPRIOR Log prior density on the Minnesota hyperparameters.
            %   lp = spec.logHyperprior(priorcoef)        % lambdas only
            %   lp = spec.logHyperprior(priorcoef, ppsi)  % lambdas + psi
            %
            %   Returns the POSITIVE log density log p(lambda1, lambda4,
            %   lambda5, [psi]). Needs only the spec's scalar hyperparameter
            %   values - no build step, no data - which is why it lives here
            %   rather than on minnesotabvarm: evaluating it should not require
            %   materialising a model.
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

            lp = 0;

            if isfield(priorcoef, "lambda1")
                lp = lp + minnesotaSpec.gammaLogPdf( ...
                    spec.lambda1, priorcoef.lambda1.k, priorcoef.lambda1.theta);
            end

            if isfield(priorcoef, "lambda4") && isfinite(spec.lambda4)
                lp = lp + minnesotaSpec.gammaLogPdf( ...
                    spec.lambda4, priorcoef.lambda4.k, priorcoef.lambda4.theta);
            end

            if isfield(priorcoef, "lambda5") && isfinite(spec.lambda5)
                lp = lp + minnesotaSpec.gammaLogPdf( ...
                    spec.lambda5, priorcoef.lambda5.k, priorcoef.lambda5.theta);
            end

            if isfield(priorcoef, "psi") && ~isempty(ppsi)
                lp = lp + sum(minnesotaSpec.invGammaLogPdf( ...
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