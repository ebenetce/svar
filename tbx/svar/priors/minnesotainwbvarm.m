classdef (Hidden) minnesotainwbvarm < semiconjugatebvarm & svar.minnesotabvarmBase & matlab.mixin.CustomDisplay
    %MINNESOTAINWBVARM Independent Normal-Wishart Minnesota prior for a BVAR.
    %
    %   A subclass of SEMICONJUGATEBVARM implementing the Independent
    %   Normal-Wishart (INW) Minnesota prior of Kadiyala & Karlsson (1997,
    %   Journal of Applied Econometrics) - the prior family used whenever a
    %   GENUINELY free own-vs-cross-variable shrinkage hyperparameter
    %   (lambda2) is required, which the conjugate Normal-Wishart prior
    %   (MINNESOTABVARM, subclassing CONJUGATEBVARM) cannot represent exactly.
    %
    %   Why this needs a different superclass, not just different moments
    %   -------------------------------------------------------------------
    %   Under the conjugate prior, Var(vec B) = Omega (x) V - a Kronecker
    %   product. The variance of coefficient (regressor r, target i) is then
    %   Omega_ii * V_rr: a product of a function of i alone and a function of
    %   r alone. "Own vs cross" is the comparison i = k vs i != k (k = the
    %   regressor's source variable) - a relationship BETWEEN the two
    %   indices, which no product of separately-indexed factors can encode,
    %   for ANY choice of V, diagonal or dense. This is not a limitation of
    %   any particular V; it is a structural ceiling of the Kronecker form
    %   itself (see MINNESOTABVARM's header for the full derivation).
    %
    %   The INW prior removes the Kronecker link entirely: V is a full
    %   (m*NumSeries)-by-(m*NumSeries) matrix (here, diagonal) indexed
    %   directly over (regressor, target) PAIRS, so the own/cross distinction
    %   can be written straight into the diagonal. The price is that the
    %   coefficient prior (Normal) and the innovations-covariance prior
    %   (Inverse-Wishart on Sigma) are no longer conjugate to the joint
    %   likelihood - hence SEMICONJUGATEBVARM, and posterior draws via Gibbs
    %   sampling (estimate returns an EMPIRICALBVARM), not a closed form.
    %
    %   Design contract (same spirit as MINNESOTABVARM, adapted to INW)
    %   ------------------------------------------------------------------
    %   * The constructor is data-free: it takes precomputed residual
    %     variances, not Y (see ESTIMATERESIDUALVARIANCES), so the object is
    %     a light, reusable recipe.
    %   * The object remembers its hyperparameters as read-only properties.
    %   * estimate() is a THIN passthrough to SEMICONJUGATEBVARM's Gibbs
    %     sampler. Set rng() before calling - the sampler is stochastic, and
    %     nothing in this class can substitute for that.
    %
    %   NOT supported (documented gaps, not silent approximations)
    %   -------------------------------------------------------------------
    %   * lambda4 / lambda5 (sum-of-coefficients, dummy-initial-observation
    %     dummies) are NOT implemented here. Folding these into MINNESOTABVARM
    %     works via a closed-form conjugate NIW update that depends on the
    %     Omega/V Kronecker link this class deliberately does not have.
    %     Threading dummy rows through SEMICONJUGATEBVARM's Gibbs sampler
    %     would require its exact internal design-matrix override mechanics,
    %     which cannot be verified without a live MATLAB session - rather
    %     than guess, this is left unimplemented.
    %   * logMarginalLikelihood is NOT provided. Gibbs draws have no closed-
    %     form marginal likelihood; obtaining one requires a separate,
    %     substantial computation (e.g. Chib's (1995) method). Hyperparameter
    %     selection for THIS class should use an out-of-sample criterion
    %     (rolling-window forecast RMSE/log score over a lambda2 grid, say),
    %     not GLPOPTIMIZEMINNESOTA, which requires the analytic marginal
    %     likelihood that only MINNESOTABVARM provides.
    %
    %   Hyperparameter mapping (own/cross target variance)
    %   -------------------------------------------------------------------
    %       Var(B_{lag l, source k -> target i}) =
    %           lambda1^2 / l^(2*lambda3) * (psi_i/psi_k)                 if i = k (own)
    %           lambda1^2 * lambda2^2 / l^(2*lambda3) * (psi_i/psi_k)     if i ~= k (cross)
    %   the standard Litterman / Kadiyala-Karlsson target (see e.g. Canova
    %   2007, Ch.10). lambda2 in (0,1] is the free cross-variable tightness
    %   this class exists to provide; lambda2 = 1 recovers the same target
    %   MINNESOTABVARM approximates via its Kronecker form (own = cross
    %   shape, only the psi ratio differs) - so lambda2 = 1 is a useful
    %   sanity check when comparing the two classes on the same data.

    properties (SetAccess = private)
        lambda2   (1,1) double   % cross-variable relative tightness (FREE here)
    end

    methods

        function obj = minnesotainwbvarm(numseries, numlags, nvp, nvp2)
            arguments
                numseries (1,1) double {mustBeInteger, mustBePositive}
                numlags   (1,1) double {mustBeInteger, mustBePositive}
                nvp.ResidualVariances (1,:) double {mustBePositive} = []
                nvp.lambda1   (1,1) double {mustBePositive}    = 0.2
                nvp.lambda2   (1,1) double {mustBePositive}    = 0.5
                nvp.lambda3   (1,1) double {mustBeNonnegative} = 1
                nvp.Vc        (1,1) double {mustBePositive}    = 1e4
                nvp.PriorMean (1,:) double = []
                nvp2.Description
                nvp2.IncludeConstant
                nvp2.IncludeTrend
                nvp2.NumPredictors
                nvp2.SeriesNames
            end

            args = namedargs2cell(nvp2);
            obj  = obj@semiconjugatebvarm(numseries, numlags, args{:});

            [residualVariances, priorMean] = obj.validateMinnesotaInputs( ...
                nvp.ResidualVariances, nvp.PriorMean, "minnesotainwbvarm");

            obj.ResidualVariances = residualVariances;
            obj.lambda1           = nvp.lambda1;
            obj.lambda2           = nvp.lambda2;
            obj.lambda3           = nvp.lambda3;
            obj.Vc                = nvp.Vc;
            obj.PriorMean         = priorMean;

            [obj.Mu, obj.V, obj.Omega, obj.DoF] = obj.buildIndependentPrior();
        end

        % function estimate
        % end

    end

    % ---- prior construction (data-free) -----------------------------------
    methods (Access = private)

        function [Mu, V, Omega, DoF] = buildIndependentPrior(obj)
            %BUILDINDEPENDENTPRIOR INW Minnesota moments from residual variances.
            %   V is the FULL (m*n)-by-(m*n) diagonal coefficient covariance
            %   (no Kronecker link to Omega), indexed directly over
            %   (regressor, target) pairs so own/cross shrinkage can differ.
            %   Coefficient layout matches conjugatebvarm/semiconjugatebvarm:
            %     vec([Phi1 ... PhiP  c  delta  B]'), an m-by-n matrix.
            n     = obj.NumSeries;
            DoF   = n + 2;              % independent IW prior on Sigma, same
            Omega = diag(obj.ResidualVariances); % minimal-informative convention as
            % MINNESOTABVARM (see its header).

            Mu = obj.buildMinnesotaPriorMean();
            V  = obj.buildIndependentCoefficientCovariance(obj.lambda2);
        end

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

    % ---- display -----------------------------------------------------------
    methods (Access = protected)

        function displayScalarObject(obj)
            disp(matlab.mixin.CustomDisplay.getSimpleHeader(obj));
            base  = {'NumSeries','P','ResidualVariances','lambda1','lambda2','lambda3', ...
                'Vc','PriorMean'};
            group = obj.minnesotaPropertyGroup(base);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj, group);
        end

    end

end
