classdef (Hidden) minnesotamniwbvarm < conjugatebvarm & svar.minnesotabvarmBase & matlab.mixin.CustomDisplay
    %MINNESOTAMNIWBVARM Conjugate Minnesota prior for a Bayesian VAR.
    %
    %   A thin subclass of CONJUGATEBVARM whose constructor materialises a
    %   Litterman-shaped Matrix-Normal-Inverse-Wishart prior from a residual
    %   variance vector (ResidualVariances) and a small set of hyperparameters.
    %   It is a *base prior* only: it never stores the estimation sample Y.
    %
    %   Design contract
    %   ---------------
    %   * The constructor is data-free. It takes precomputed residual
    %     variances, not Y, so the object is a light, reusable recipe and the
    %     GLP-style optimiser can rebuild it cheaply for each candidate
    %     hyperparameter set without refitting anything.
    %   * The object remembers its hyperparameters as read-only properties,
    %     so it can be inspected, displayed, and re-materialised with one
    %     field changed.
    %   * Sum-of-coefficients (lambda4) and dummy-initial-observation
    %     (lambda5) priors are data-dependent. They are NOT baked into the
    %     constructor. They are applied - via the same conjugate NIW update
    %     used for the posterior - at the one place data legitimately enters:
    %     estimate / logMarginalLikelihood / simulate / forecast. A single
    %     private helper (applyDummies) builds the dummy-augmented prior so
    %     estimation and the marginal likelihood always use an identical
    %     augmented prior.
    %   * lambda4 = lambda5 = Inf means "dummy off". The switch is isfinite,
    %     not a zero pseudo-observation, so an Inf lambda contributes no rows
    %     at all (a zero row would still inflate the degrees of freedom).
    %
    %   Hyperparameter mapping
    %   ----------------------
    %       lambda1  overall Minnesota tightness           (paper lambda)
    %       lambda3  lag-decay exponent; variance ~ 1/l^(2*lambda3), paper = 1
    %       lambda4  sum-of-coefficients tightness         (paper mu; Inf = off)
    %       lambda5  dummy-initial-observation tightness   (paper delta; Inf = off)
    %       Vc       prior variance of the constant / trend
    %       ResidualVariances  per-series residual variances (the prior scale)
    %
    %   Note on the fixed lambda2. The classic Litterman prior has a second
    %   hyperparameter (cross-variable relative tightness). It is deliberately
    %   absent here: the conjugate Kronecker structure V = KK (x) inv(Psi)
    %   cannot represent a separate own-vs-cross tightness, so lambda2 is
    %   implicitly 1. Representing it would break conjugacy - and conjugacy is
    %   exactly what makes the analytic marginal likelihood (and the GLP
    %   tuner) possible. If you need a free lambda2, that belongs in a
    %   semiconjugate variant, not this class.

    properties (SetAccess = private)
        lambda2   (1,1) double = 1 % cross-variable relative tightness (pinned to 1)
        lambda4   (1,1) double     % sum-of-coefficients tightness (Inf = off)
        lambda5   (1,1) double     % dummy-initial-observation tightness (Inf = off)
    end

    methods

        function obj = minnesotamniwbvarm(numseries, numlags, residualVariances, nvp, nvp2)
            arguments
                numseries (1,1) double {mustBeInteger, mustBePositive}
                numlags   (1,1) double {mustBeInteger, mustBePositive}
                residualVariances (1,:) double {mustBePositive}
                nvp.lambda1   (1,1) double {mustBePositive}    = 0.2
                nvp.lambda3   (1,1) double {mustBeNonnegative} = 1
                nvp.lambda4   (1,1) double {mustBePositive}    = Inf
                nvp.lambda5   (1,1) double {mustBePositive}    = Inf
                nvp.Vc        (1,1) double {mustBePositive}    = 1e4
                nvp.PriorMean (1,:) double = ones(1, numseries);
                nvp2.Description
                nvp2.IncludeConstant
                nvp2.IncludeTrend
                nvp2.NumPredictors
                nvp2.SeriesNames
            end

            % Delegate the structural set-up (SeriesNames, exogenous layout,
            % NumSeries/P bookkeeping) to the conjugate superclass.
            args = namedargs2cell(nvp2);
            obj  = obj@conjugatebvarm(numseries, numlags, args{:});

            priorMean = obj.validateMinnesotaInputs( ...
                residualVariances, nvp.PriorMean, "minnesotamniwbvarm");

            % Store the hyperparameters (read-only from here on).
            obj.ResidualVariances = residualVariances;
            obj.lambda1           = nvp.lambda1;
            obj.lambda3           = nvp.lambda3;
            obj.lambda4           = nvp.lambda4;
            obj.lambda5           = nvp.lambda5;
            obj.Vc                = nvp.Vc;
            obj.PriorMean         = priorMean;

            % Materialise the BASE prior moments (no dummies, no data).
            [obj.Mu, obj.V, obj.Omega, obj.DoF] = obj.buildBasePrior();
        end

    end
    % ---- public API ------------------------------------------------------
    methods

        function [Posterior, Summary] = estimate(obj, Y, opts)
            %ESTIMATE Analytic conjugate posterior, dummy observations included.
            arguments
                obj
                Y double {mustBeNonempty}
                opts.Display
                opts.X
                opts.Y0
            end

            prior = obj.applyDummies(Y);            % dummy-augmented prior
            args  = namedargs2cell(opts);
            [MN, Summary] = estimate@conjugatebvarm(prior, Y, args{:});

            % Return a plain conjugatebvarm posterior (it is no longer a
            % Minnesota prior, so do not pretend it is one).
            Posterior = conjugatebvarm(MN.NumSeries, MN.P, ...
                Description     = MN.Description, ...
                SeriesNames     = MN.SeriesNames, ...
                IncludeConstant = MN.IncludeConstant, ...
                IncludeTrend    = MN.IncludeTrend, ...
                NumPredictors   = MN.NumPredictors, ...
                Mu = MN.Mu, V = MN.V, Omega = MN.Omega, DoF = MN.DoF);
        end

        function varargout = simulate(obj, Y, opts)
            %SIMULATE Draw from the (dummy-augmented) prior given data.
            arguments
                obj
                Y double {mustBeNonempty}
                opts.NumDraws
                opts.X
                opts.Y0
            end
            prior = obj.applyDummies(Y);
            args  = namedargs2cell(opts);
            [varargout{1:nargout}] = simulate@conjugatebvarm(prior, Y, args{:});
        end

        function varargout = forecast(obj, Y, numperiods)
            %FORECAST Forecast from the (dummy-augmented) prior given data.
            arguments
                obj
                Y          double {mustBeNonempty}
                numperiods (1,1) double {mustBeInteger, mustBePositive}
            end
            prior = obj.applyDummies(Y);
            [varargout{1:nargout}] = forecast@conjugatebvarm(prior, numperiods, Y);
        end

        function varargout = simsmooth(obj, Y, varargin)
            %SIMSMOOTH Simulation smoother on the (dummy-augmented) prior.
            arguments
                obj
                Y double {mustBeNonempty}
            end
            arguments (Repeating)
                varargin
            end
            prior = obj.applyDummies(Y);
            [varargout{1:nargout}] = simsmooth@conjugatebvarm(prior, Y, varargin{:});
        end

    end

    % ---- prior construction (data-free) ---------------------------------
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

    end

    % ---- data-dependent pieces (dummies, regression matrices) -----------
    methods

        function prior = applyDummies(obj, Y)
            %APPLYDUMMIES Fold any active dummy priors into the base prior.
            %   Returns a copy of OBJ with Mu/V/Omega/DoF updated via the
            %   conjugate NIW update. Shared by estimate and the marginal
            %   likelihood so both always see the identical augmented prior.
            prior = obj;
            [XDummy, YDummy] = obj.dummyMatrices(Y);
            if isempty(YDummy)
                return
            end
            [prior.Mu, prior.V, prior.Omega, prior.DoF] = ...
                minnesotamniwbvarm.updateNIW(obj.Mu, obj.V, obj.Omega, obj.DoF, ...
                                         XDummy, YDummy, obj.NumSeries);
        end

        function [XDummy, YDummy] = dummyMatrices(obj, Y)
            %DUMMYMATRICES Sum-of-coefficients and initial-observation dummies.
            %   Built only for finite lambda4 / lambda5 (Inf = off). ybar is
            %   the mean of the first P observations of each series.
            obj.validateMinnesotaData(Y, "minnesotamniwbvarm");
            n  = obj.NumSeries;
            P  = obj.P;
            mm = obj.m;

            XDummy = zeros(0, mm);
            YDummy = zeros(0, n);
            if ~isfinite(obj.lambda4) && ~isfinite(obj.lambda5)
                return
            end

            ybar = mean(Y(1:P, :), 1);                 % 1-by-n

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

    end

    % ---- shared conjugate NIW algebra -----------------------------------
    methods (Static, Access = private)

        function [MuOut, VOut, OmegaOut, DoFOut] = updateNIW(Mu, V, Omega, DoF, X, Yobs, n)
            %UPDATENIW One conjugate Normal-Inverse-Wishart update.
            %   Folds the pseudo-/real observations (X, Yobs) into the prior.
            %   The same update produces the posterior and the dummy-augmented
            %   prior; factoring it here keeps a single, tested code path.
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

            base  = {'NumSeries','P','ResidualVariances','lambda1','lambda3','lambda4','lambda5','Vc','PriorMean'};
            
            group = obj.minnesotaPropertyGroup(base);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj, group);
        end

    end

end
