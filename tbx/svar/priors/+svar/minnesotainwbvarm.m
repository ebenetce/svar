classdef minnesotainwbvarm < semiconjugatebvarm & svar.minnesotabvarmBase & matlab.mixin.CustomDisplay
    %MINNESOTAINWBVARM Independent Normal-Wishart Minnesota prior for a BVAR.
    %
    %   PriorMdl = SVAR.MINNESOTAINWBVARM(NUMSERIES,NUMLAGS,Y) creates the
    %   Independent Normal-Wishart (INW) Minnesota prior of Kadiyala &
    %   Karlsson (1997, Journal of Applied Econometrics) for the sample Y.
    %   Prefer the MINNESOTABVARM front door, which reaches this class as
    %   Type="inw". This is the family to use whenever a GENUINELY free
    %   own-vs-cross-variable shrinkage hyperparameter (lambda2) is required,
    %   which the conjugate prior cannot represent exactly.
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
    %   any particular V; it is a structural ceiling of the Kronecker form.
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
    %   Set RNG before estimating or simulating: the sampler is stochastic,
    %   and nothing in this class can substitute for that.
    %
    %   NOT supported (documented gaps, not silent approximations)
    %   -------------------------------------------------------------------
    %   * lambda4 / lambda5 (sum-of-coefficients, dummy-initial-observation
    %     dummies) are NOT implemented here. Folding these into the conjugate
    %     family works via a closed-form NIW update that depends on the
    %     Omega/V Kronecker link this class deliberately does not have.
    %     Threading dummy rows through SEMICONJUGATEBVARM's Gibbs sampler
    %     would require its exact internal design-matrix mechanics; rather
    %     than guess, this is left unimplemented.
    %   * logMarginalLikelihood is NOT provided. Gibbs draws have no closed-
    %     form marginal likelihood; obtaining one requires a separate,
    %     substantial computation (e.g. Chib's (1995) method). Hyperparameter
    %     selection for THIS class should use an out-of-sample criterion
    %     (rolling-window forecast RMSE / log score over a lambda2 grid, say),
    %     not GLP, which requires the analytic marginal likelihood only the
    %     conjugate family provides.
    %
    %   Hyperparameter mapping (own/cross target variance)
    %   -------------------------------------------------------------------
    %       Var(B_{lag l, source k -> target i}) =
    %           lambda1^2 / l^(2*lambda3) * (psi_i/psi_k)                 if i = k (own)
    %           lambda1^2 * lambda2^2 / l^(2*lambda3) * (psi_i/psi_k)     if i ~= k (cross)
    %   the standard Litterman / Kadiyala-Karlsson target (see e.g. Canova
    %   2007, Ch.10). lambda2 in (0,1] is the free cross-variable tightness
    %   this class exists to provide; lambda2 = 1 recovers the same target
    %   the conjugate family approximates via its Kronecker form (own = cross
    %   shape, only the psi ratio differs) - a useful sanity check when
    %   comparing the two on the same data.
    %
    %   See also MINNESOTABVARM, SVAR.MINNESOTAMNIWBVARM, SEMICONJUGATEBVARM.

    properties (SetAccess = private)
        lambda2   (1,1) double   % cross-variable relative tightness (FREE here)
    end

    methods

        function obj = minnesotainwbvarm(numseries, numlags, Y, nvp, nvp2)
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

            [Y, nvp2] = svar.minnesotainwbvarm.resolveSample( ...
                Y, numseries, numlags, nvp2, "minnesotainwbvarm");

            args = namedargs2cell(nvp2);
            obj  = obj@semiconjugatebvarm(numseries, numlags, args{:});

            obj = obj.configureMinnesota(Y, numlags, nvp, "minnesotainwbvarm");
            obj.lambda2 = nvp.lambda2;

            [obj.Mu, obj.V, obj.Omega, obj.DoF] = obj.buildIndependentPrior();
        end

    end

    % ---- public API ------------------------------------------------------
    % Thin overrides that default the sample to the stored one. estimate here
    % is a passthrough to the inherited Gibbs sampler.
    methods

        function varargout = estimate(obj, varargin)
            %ESTIMATE Gibbs posterior for the stored sample (EMPIRICALBVARM).
            [varargout{1:nargout}] = estimate@semiconjugatebvarm(obj, obj.Y, varargin{:});
        end

        function varargout = simulate(obj, varargin)
            %SIMULATE Draw coefficients and covariance given the stored sample.
            [varargout{1:nargout}] = simulate@semiconjugatebvarm(obj, obj.Y, varargin{:});
        end

        function varargout = forecast(obj, numperiods)
            %FORECAST Forecast responses beyond the stored sample.            
            [varargout{1:nargout}] = ...
                forecast@semiconjugatebvarm(obj, numperiods, obj.Y);
        end

        function varargout = simsmooth(obj, varargin)
            %SIMSMOOTH Simulation smoother, defaulting to the stored sample.
            %   Not sample-checked - see SVAR.MINNESOTAMNIWBVARM/SIMSMOOTH for
            %   why the inherited machinery must be allowed a derived sample.
            [varargout{1:nargout}] = simsmooth@semiconjugatebvarm(obj, obj.Y, varargin{:});
        end

    end

    % ---- prior construction ----------------------------------------------
    methods (Access = private)

        function [Mu, V, Omega, DoF] = buildIndependentPrior(obj)
            %BUILDINDEPENDENTPRIOR INW Minnesota moments from residual variances.
            %   Coefficient layout matches conjugatebvarm/semiconjugatebvarm:
            %     vec([Phi1 ... PhiP  c  delta  B]'), an m-by-n matrix.
            n     = obj.NumSeries;
            DoF   = n + 2;                        % independent IW prior on Sigma,
            Omega = diag(obj.ResidualVariances);  % same minimal-informative
                                                  % convention as the conjugate family.

            Mu = obj.buildMinnesotaPriorMean();
            V  = obj.buildIndependentCoefficientCovariance(obj.lambda2);
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
