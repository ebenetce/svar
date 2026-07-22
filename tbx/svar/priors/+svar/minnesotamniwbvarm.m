classdef minnesotamniwbvarm < conjugatebvarm & svar.minnesotabvarmBase & matlab.mixin.CustomDisplay
    %MINNESOTAMNIWBVARM Conjugate Minnesota prior for a Bayesian VAR.
    %
    %   PriorMdl = SVAR.MINNESOTAMNIWBVARM(NUMSERIES,NUMLAGS,Y) creates a
    %   Litterman-shaped Matrix-Normal-Inverse-Wishart prior for the sample Y.
    %   Prefer the MINNESOTABVARM front door, which reaches this class as
    %   Type="mniw" (the default).
    %
    %   Everything data-dependent is resolved in the constructor: the
    %   residual-variance scale Psi (see SVAR.MINNESOTABVARMBASE) and the
    %   sum-of-coefficients / dummy-initial-observation priors. The object
    %   that comes out is the FULL prior - Mu/V/Omega/DoF already carry the
    %   dummy observations.
    %
    %   That is the point of taking Y up front. A data-free object that folds
    %   dummies in later forces every consumer (estimate, simulate, forecast,
    %   simsmooth, logMarginalLikelihood) to remember to apply them, and a
    %   consumer that forgets silently drops lambda4/lambda5 rather than
    %   failing. Building the augmented prior once makes that unreachable.
    %
    %   Hyperparameter mapping
    %   ----------------------
    %       lambda1  overall Minnesota tightness           (paper lambda)
    %       lambda3  lag-decay exponent; variance ~ 1/l^(2*lambda3), paper = 1
    %       lambda4  sum-of-coefficients tightness         (paper mu; Inf = off)
    %       lambda5  dummy-initial-observation tightness   (paper delta; Inf = off)
    %       Vc       prior variance of the constant / trend
    %       Psi      per-series residual variances (the prior scale)
    %
    %   lambda4 = lambda5 = Inf means "dummy off". The switch is ISFINITE, not
    %   a zero pseudo-observation, so an Inf lambda contributes no rows at all
    %   (a zero row would still inflate the degrees of freedom).
    %
    %   Note on the fixed lambda2. The classic Litterman prior has a second
    %   hyperparameter (cross-variable relative tightness). It is deliberately
    %   absent here: the conjugate Kronecker structure V = KK (x) inv(Psi)
    %   cannot represent a separate own-vs-cross tightness, so lambda2 is
    %   implicitly 1. Representing it would break conjugacy - and conjugacy is
    %   exactly what makes the analytic marginal likelihood (and GLP) possible.
    %   If you need a free lambda2, use SVAR.MINNESOTAINWBVARM.
    %
    %   See also MINNESOTABVARM, GLP, LOGMARGINALLIKELIHOOD, CONJUGATEBVARM.

    properties (SetAccess = private)
        lambda2   (1,1) double = 1  % cross-variable relative tightness (pinned to 1)
        lambda4   (1,1) double      % sum-of-coefficients tightness (Inf = off)
        lambda5   (1,1) double      % dummy-initial-observation tightness (Inf = off)
    end

    properties (Dependent, SetAccess = private)
        NumDummyObservations        % pseudo-observations folded into the prior
    end

    methods

        function obj = minnesotamniwbvarm(numseries, numlags, Y, nvp, nvp2)
            arguments
                numseries (1,1) double {mustBeInteger, mustBePositive}
                numlags   (1,1) double {mustBeInteger, mustBePositive}
                Y                {mustBeNonempty}
                nvp.Psi           = "exact"
                nvp.lambda1   (1,1) double {mustBePositive}    = 0.2
                nvp.lambda3   (1,1) double {mustBeNonnegative} = 1
                nvp.lambda4   (1,1) double {mustBePositive}    = Inf
                nvp.lambda5   (1,1) double {mustBePositive}    = Inf
                nvp.Vc        (1,1) double {mustBePositive}    = 1e4
                nvp.PriorMean (1,:) double = ones(1, numseries)
                nvp2.Description
                nvp2.IncludeConstant
                nvp2.IncludeTrend
                nvp2.NumPredictors
                nvp2.SeriesNames
            end

            [Y, nvp2] = svar.minnesotamniwbvarm.resolveSample( ...
                Y, numseries, numlags, nvp2, "minnesotamniwbvarm");

            % Delegate the structural set-up (SeriesNames, exogenous layout,
            % NumSeries/P bookkeeping) to the conjugate superclass.
            args = namedargs2cell(nvp2);
            obj  = obj@conjugatebvarm(numseries, numlags, args{:});

            obj = obj.configureMinnesota(Y, numlags, nvp, "minnesotamniwbvarm");
            obj.lambda4 = nvp.lambda4;
            obj.lambda5 = nvp.lambda5;

            % Materialise the FULL prior: Litterman moments, then any active
            % dummy observations folded in through the conjugate NIW update.
            [obj.Mu, obj.V, obj.Omega, obj.DoF] = obj.buildBasePrior();
            [obj.Mu, obj.V, obj.Omega, obj.DoF] = obj.applyDummies();
        end

        function value = get.NumDummyObservations(obj)
            [~, YDummy] = obj.dummyMatrices();
            value = size(YDummy, 1);
        end

    end

    % ---- public API ------------------------------------------------------
    % These overrides exist only to default the sample to the stored one and
    % to reject a mismatched sample. The prior itself needs no late patching.
    methods

        function [Posterior, Summary] = estimate(obj, varargin)
            %ESTIMATE Analytic conjugate posterior for the stored sample.
      
            Mdl = toConjugate(obj);
            [Posterior, Summary] = Mdl.estimate(obj.Y, varargin{:});

        end

        function varargout = simulate(obj, varargin)
            %SIMULATE Draw coefficients and covariance given the stored sample.

            [varargout{1:nargout}] = simulate@conjugatebvarm(obj, obj.Y, varargin{:});
        end

        function varargout = forecast(obj, numperiods)
            %FORECAST Forecast responses beyond the stored sample.
            arguments
                obj        (1,1) svar.minnesotamniwbvarm
                numperiods (1,1) double {mustBeInteger, mustBePositive}
            end
            Mdl = toConjugate(obj);
            [varargout{1:nargout}] = Mdl.forecast(numperiods, obj.Y);
        end

        function varargout = simsmooth(obj, varargin)
            %SIMSMOOTH Simulation smoother, defaulting to the stored sample.
            %   Deliberately NOT sample-checked. simsmooth is the primitive
            %   the inherited machinery calls internally - bvar/forecast in
            %   particular passes the stored sample extended with NaN rows
            %   over the forecast horizon - so a strict equality check here
            %   would reject the toolbox's own legitimate calls. The check
            %   belongs on the user-facing entry points above.
            [varargout{1:nargout}] = simsmooth@conjugatebvarm(obj, obj.Y, varargin{:});
        end

    end

    % ---- prior construction ---------------------------------------------
    methods (Access = private)

        function [Mu, V, Omega, DoF] = buildBasePrior(obj)
            %BUILDBASEPRIOR Litterman moments from residual variances + hyperparameters.
            %   Coefficient layout matches conjugatebvarm:
            %     vec([Phi1 ... PhiP  c  delta  B]'), i.e. an m-by-n matrix,
            %     lag blocks first, then constant, then trend, then predictors.
            n   = obj.NumSeries;
            P   = obj.P;
            mm  = obj.m;
            psi = obj.ResidualVariances(:);

            % DoF fixed at n+2: the minimal proper IW with a defined mean, the
            % standard Minnesota / GLP choice. Overrides the inherited n+10
            % conjugatebvarm default, which V (below) is NOT built for.
            DoF = n + 2;

            Omega = diag(psi);

            % Prior mean: own first lag only.
            Mu = obj.buildMinnesotaPriorMean();

            % Prior covariance factor V (diagonal).
            %   With DoF = n+2 the (DoF-n-1) scaling is 1 by construction, so
            %   the marginal prior coefficient variance equals the Minnesota
            %   target directly: Var(vec B) = Omega/(DoF-n-1) (x) V = Omega (x) V.
            %   Lag block:  V(lag l, regressor j) = lambda1^2 / (l^(2*l3) * psi_j).
            lagScale = obj.lambda1^2 ./ ((1:P)'.^(2*obj.lambda3));   % P-by-1
            KK       = diag(lagScale);                 % P-by-P
            Vlag     = kron(KK, diag(1 ./ psi));       % (P*n)-by-(P*n)

            vDiag                 = zeros(mm, 1);
            vDiag(1:P*n)          = diag(Vlag);
            vDiag(P*n+1 : mm)     = obj.Vc;            % constant, trend, predictors
            V = diag(vDiag);
        end

        function [Mu, V, Omega, DoF] = applyDummies(obj)
            %APPLYDUMMIES Fold any active dummy priors into the base prior.
            [Mu, V, Omega, DoF] = deal(obj.Mu, obj.V, obj.Omega, obj.DoF);
            [XDummy, YDummy] = obj.dummyMatrices();
            if isempty(YDummy)
                return
            end
            [Mu, V, Omega, DoF] = svar.minnesotamniwbvarm.updateNIW( ...
                Mu, V, Omega, DoF, XDummy, YDummy, obj.NumSeries);
        end

        function [XDummy, YDummy] = dummyMatrices(obj)
            %DUMMYMATRICES Sum-of-coefficients and initial-observation dummies.
            %   Built only for finite lambda4 / lambda5 (Inf = off). ybar is
            %   the mean of the first P observations of each series.
            n  = obj.NumSeries;
            P  = obj.P;
            mm = obj.m;

            XDummy = zeros(0, mm);
            YDummy = zeros(0, n);
            if ~isfinite(obj.lambda4) && ~isfinite(obj.lambda5)
                return
            end

            ybar = mean(obj.Y(1:P, :), 1);             % 1-by-n

            % Sum-of-coefficients: n rows, one per series.
            if isfinite(obj.lambda4)
                socY  = diag(ybar) / obj.lambda4;      % n-by-n
                Xsoc  = zeros(n, mm);
                for lag = 1:P
                    cols = (lag-1)*n + (1:n);
                    Xsoc(:, cols) = socY;
                end
                XDummy = [XDummy; Xsoc];
                YDummy = [YDummy; socY];
            end

            % Dummy initial observation (co-persistence): a single row.
            if isfinite(obj.lambda5)
                dioY  = ybar / obj.lambda5;            % 1-by-n
                Xdio  = zeros(1, mm);
                for lag = 1:P
                    cols = (lag-1)*n + (1:n);
                    Xdio(:, cols) = dioY;
                end
                if obj.IncludeConstant
                    Xdio(1, P*n + 1) = 1 / obj.lambda5;
                end
                % Trend / predictor columns left at 0 (no co-persistence dummy).
                XDummy = [XDummy; Xdio];
                YDummy = [YDummy; dioY];
            end
        end

        function Mdl = toConjugate(obj)

            Mdl = conjugatebvarm(obj.NumSeries, obj.P, ...
                Description     = obj.Description, ...
                SeriesNames     = obj.SeriesNames, ...
                IncludeConstant = obj.IncludeConstant, ...
                IncludeTrend    = obj.IncludeTrend, ...
                NumPredictors   = obj.NumPredictors, ...
                Mu = obj.Mu, V = obj.V, Omega = obj.Omega, DoF = obj.DoF);

        end

    end

    % ---- shared conjugate NIW algebra -----------------------------------
    methods (Static, Access = private)

        function [MuOut, VOut, OmegaOut, DoFOut] = updateNIW(Mu, V, Omega, DoF, X, Yobs, n)
            %UPDATENIW One conjugate Normal-Inverse-Wishart update.
            %   Folds the pseudo-observations (X, Yobs) into the prior.
            k      = size(V, 1);
            B0     = reshape(Mu, k, n);
            V0inv  = V \ eye(k);
            prec   = V0inv + X'*X;
            B      = prec \ (V0inv*B0 + X'*Yobs);
            Vnew   = prec \ eye(k);
            OmegaN = Omega + Yobs'*Yobs + B0'*V0inv*B0 - B'*prec*B;

            MuOut    = B(:);
            VOut     = (Vnew + Vnew')/2;
            OmegaOut = (OmegaN + OmegaN')/2;
            DoFOut   = DoF + size(Yobs, 1);
        end

    end

    % ---- display ---------------------------------------------------------
    methods (Access = protected)

        function displayScalarObject(obj)
            disp(matlab.mixin.CustomDisplay.getSimpleHeader(obj));

            base = {'NumSeries','P','ResidualVariances','lambda1','lambda3', ...
                'lambda4','lambda5','Vc','PriorMean','NumDummyObservations'};

            group = obj.minnesotaPropertyGroup(base);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj, group);
        end

    end

end