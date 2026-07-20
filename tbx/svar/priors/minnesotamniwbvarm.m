classdef (Hidden) minnesotamniwbvarm < conjugatebvarm & minnesotabvarm & matlab.mixin.CustomDisplay
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

        function obj = minnesotamniwbvarm(numseries, numlags, nvp, nvp2)
            arguments
                numseries (1,1) double {mustBeInteger, mustBePositive}
                numlags   (1,1) double {mustBeInteger, mustBePositive}
                nvp.ResidualVariances (1,:) double {mustBePositive} = []
                nvp.lambda1   (1,1) double {mustBePositive}    = 0.2
                nvp.lambda3   (1,1) double {mustBeNonnegative} = 1
                nvp.lambda4   (1,1) double {mustBePositive}    = Inf
                nvp.lambda5   (1,1) double {mustBePositive}    = Inf
                nvp.Vc        (1,1) double {mustBePositive}    = 1e4
                nvp.PriorMean (1,:) double = []
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

            [residualVariances, priorMean] = obj.validateMinnesotaInputs( ...
                nvp.ResidualVariances, nvp.PriorMean, "minnesotamniwbvarm");

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

        function [logML, details] = logMarginalLikelihood(obj, Y)
            %LOGMARGINALLIKELIHOOD  log p(Y | hyperparameters).
            %   With dummy observations this is evaluated against the
            %   dummy-augmented prior using the REAL data only, which equals
            %   p(Y, dummies) / p(dummies): the dummies' own evidence cancels.
            %   This is the correct objective for hyperparameter selection;
            %   feeding the dummy rows in as extra data would instead return
            %   p(Y, dummies) and bias the tuner. This returns the PURE
            %   marginal likelihood - any hyperprior on the lambdas belongs in
            %   the optimiser layer (see logHyperprior), not here.
            arguments
                obj
                Y double {mustBeNonempty}
            end
            prior          = obj.applyDummies(Y);
            [X, YResponse] = obj.regressionMatrices(Y);
            [logML, details] = minnesotamniwbvarm.conjugateLogML( ...
                prior.Mu, prior.V, prior.Omega, prior.DoF, X, YResponse, obj.NumSeries);
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
    methods (Access = private)

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
            obj.validateData(Y);
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

        function [X, YResponse] = regressionMatrices(obj, Y)
            %REGRESSIONMATRICES VAR design matrix consistent with the prior layout.
            obj.validateData(Y);
            n      = obj.NumSeries;
            P      = obj.P;
            mm     = obj.m;
            numObs = size(Y, 1);
            T      = numObs - P;

            X = zeros(T, mm);
            for lag = 1:P
                cols = (lag-1)*n + (1:n);
                X(:, cols) = Y((P + 1 - lag):(numObs - lag), :);
            end
            col = P*n;
            if obj.IncludeConstant
                col = col + 1;
                X(:, col) = 1;
            end
            if obj.IncludeTrend
                col = col + 1;
                X(:, col) = (1:T)';
            end
            if obj.NumPredictors > 0
                error("minnesotamniwbvarm:predictorsUnsupported", ...
                    ["The marginal likelihood / dummy path needs the exogenous " ...
                     "regressors, which this prior does not store. Estimate with " ...
                     "the X name-value argument instead."]);
            end

            YResponse = Y((P + 1):numObs, :);
        end

        function validateData(obj, Y)
            if size(Y, 2) ~= obj.NumSeries
                error("minnesotamniwbvarm:invalidData", ...
                    "Y must have %d columns, one per series.", obj.NumSeries);
            end
            if size(Y, 1) <= obj.P
                error("minnesotamniwbvarm:invalidData", ...
                    "Y must have more rows than the lag order P = %d.", obj.P);
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

        function [logML, details] = conjugateLogML(Mu, V, Omega, DoF, X, Yobs, n)
            %CONJUGATELOGML Analytic log marginal likelihood for the conjugate VAR.
            %   Evaluated against the supplied (dummy-augmented) prior using the
            %   real data (X, Yobs) only.
            %
            %   Numerics: uses the stable ratio-of-determinants ("eig + 1")
            %   form. The two log-det differences that appear in the marginal
            %   likelihood are computed as log det(I + M) with M >= 0, so they
            %   stay finite even when the prior V or Omega is near-singular -
            %   which a hyperparameter optimiser will inevitably probe. The
            %   whole computation is guarded: any breakdown (e.g. chol of a
            %   collapsed prior) returns logML = -Inf, which an optimiser reads
            %   as "this hyperparameter set is very bad" rather than throwing.
            k        = size(V, 1);
            numObs   = size(Yobs, 1);
            priorDoF = DoF;
            postDoF  = DoF + numObs;

            details = struct();

            try
                B0    = reshape(Mu, k, n);
                V0inv = V \ eye(k);
                prec  = V0inv + X'*X;
                Bn    = prec \ (V0inv*B0 + X'*Yobs);

                % Posterior scale as prior scale plus a manifestly PSD
                % increment (residual SS at the posterior mean + prior-mean
                % shift). Avoids the cancellation-prone
                %   Omega + Y'Y + B0'V0inv B0 - Bn'prec Bn.
                resid = Yobs - X*Bn;
                shift = Bn - B0;
                incr  = resid'*resid + shift'*(V0inv*shift);
                incr  = (incr + incr')/2;                       % PSD

                % logdet(OmegaN) - logdet(Omega) = sum log( eig(inv(Omega)*incr) + 1 ).
                % Lo Lo' = Omega  =>  bbb = Lo^{-1} incr Lo^{-T} is symmetric
                % with eig(bbb) = eig(inv(Omega) incr) >= 0.
                Lo   = chol((Omega + Omega')/2, 'lower');
                bbb  = Lo \ incr / Lo';
                eb   = real(eig((bbb + bbb')/2));
                eb(eb < 0) = 0;
                sumOmegaRatio = sum(log(eb + 1));
                logDetOmega0  = 2*sum(log(diag(Lo)));           % reuse the factor

                % logdet(Vn) - logdet(V0) = -sum log( eig(X'X * V0) + 1 ).
                % D D' = V0  =>  aaa = D'(X'X)D is symmetric, eig = eig(X'X V0).
                D    = chol((V + V')/2, 'lower');
                aaa  = D'*(X'*X)*D;
                ea   = real(eig((aaa + aaa')/2));
                ea(ea < 0) = 0;
                sumVRatio = sum(log(ea + 1));

                logML = -0.5*numObs*n*log(pi) ...
                    + minnesotamniwbvarm.logMvGamma(0.5*postDoF,  n) ...
                    - minnesotamniwbvarm.logMvGamma(0.5*priorDoF, n) ...
                    - 0.5*numObs*logDetOmega0 ...
                    - 0.5*postDoF*sumOmegaRatio ...
                    - 0.5*n*sumVRatio;

                if nargout > 1
                    Vn     = prec \ eye(k);
                    OmegaN = (Omega + incr + (Omega + incr)')/2;
                    details = struct( ...
                        "NumObservations", numObs, ...
                        "PriorDoF",        priorDoF, ...
                        "PosteriorDoF",    postDoF, ...
                        "PosteriorMu",     Bn(:), ...
                        "PosteriorV",      (Vn + Vn')/2, ...
                        "PosteriorOmega",  OmegaN);
                end
            catch
                logML   = -Inf;
                details = struct("NumObservations", numObs, ...
                                 "PriorDoF", priorDoF, "PosteriorDoF", postDoF);
            end
        end

        function value = logMvGamma(a, dimension)
            %LOGMVGAMMA Log of the multivariate gamma function.
            j     = 1:dimension;
            value = dimension*(dimension - 1)*0.25*log(pi) ...
                  + sum(gammaln(a + 0.5*(1 - j)));
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
