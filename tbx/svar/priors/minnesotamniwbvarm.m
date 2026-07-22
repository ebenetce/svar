classdef minnesotamniwbvarm < conjugatebvarm & svar.minnesotabvarmBase & matlab.mixin.CustomDisplay
    %MINNESOTAMNIWBVARM Conjugate Minnesota prior for a Bayesian VAR.
    %
    %   PriorMdl = MINNESOTAMNIWBVARM(NUMSERIES,NUMLAGS,Y) creates a
    %   Litterman-shaped Matrix-Normal-Inverse-Wishart prior for the sample Y,
    %   with the residual-variance scale estimated from Y.
    %
    %   PriorMdl = MINNESOTAMNIWBVARM(...,Psi=PSI) sets the residual-variance
    %   scale. PSI is either a 1-by-NUMSERIES vector of variances or one of
    %   "exact" (default) / "conditional", naming the estimator that
    %   ESTIMATERESIDUALVARIANCES applies to Y. Pass PSI numerically inside a
    %   tuning loop: the string form refits an AR per series on every call.
    %
    %   PriorMdl = MINNESOTAMNIWBVARM(...,Name=Value) sets hyperparameters
    %   lambda1, lambda3, lambda4, lambda5, Vc, and PriorMean, plus any
    %   CONJUGATEBVARM option (SeriesNames, IncludeConstant, IncludeTrend,
    %   NumPredictors, Description).
    %
    %   Data-in-the-constructor contract
    %   --------------------------------
    %   Y is required and retained. Both data-dependent pieces of this prior -
    %   the residual-variance scale and the sum-of-coefficients /
    %   dummy-initial-observation priors - are resolved once, here, so the
    %   object that comes out is the FULL prior: Mu/V/Omega/DoF already carry
    %   the dummy observations.
    %
    %   That is the whole point of taking Y up front. The alternative - a
    %   data-free object that folds dummies in later - forces every consumer
    %   (estimate, simulate, forecast, simsmooth, logMarginalLikelihood) to
    %   remember to apply them, and a consumer that forgets silently drops
    %   lambda4/lambda5 rather than failing. Building the augmented prior once
    %   makes that class of bug unreachable.
    %
    %   Because Y is stored, ESTIMATE/SIMULATE/FORECAST/SIMSMOOTH and
    %   LOGMARGINALLIKELIHOOD may be called with no data argument. Supplying a
    %   sample that differs from the stored one is an error, not a silent
    %   re-fit: Psi and the dummy rows were derived from the stored Y, so a
    %   different sample needs a different prior object.
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
    %   If you need a free lambda2, use MINNESOTAINWBVARM.
    %
    %   See also GLP, ESTIMATERESIDUALVARIANCES, LOGMARGINALLIKELIHOOD,
    %   CONJUGATEBVARM.

    properties (SetAccess = private)
        Y         double            % estimation sample the prior was built from
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

            % Resolve the sample before the superclass call: a tabular Y also
            % supplies SeriesNames unless the caller named them explicitly.
            [Y, nvp2] = minnesotamniwbvarm.resolveSample(Y, numseries, numlags, nvp2);

            % Delegate the structural set-up (SeriesNames, exogenous layout,
            % NumSeries/P bookkeeping) to the conjugate superclass.
            args = namedargs2cell(nvp2);
            obj  = obj@conjugatebvarm(numseries, numlags, args{:});

            residualVariances = obj.resolvePsi(nvp.Psi, Y, numlags);
            priorMean = obj.validateMinnesotaInputs( ...
                residualVariances, nvp.PriorMean, "minnesotamniwbvarm");

            % Store the sample and hyperparameters (read-only from here on).
            obj.Y                 = Y;
            obj.ResidualVariances = residualVariances;
            obj.lambda1           = nvp.lambda1;
            obj.lambda3           = nvp.lambda3;
            obj.lambda4           = nvp.lambda4;
            obj.lambda5           = nvp.lambda5;
            obj.Vc                = nvp.Vc;
            obj.PriorMean         = priorMean;

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

        function [Posterior, Summary] = estimate(obj, Y, opts)
            %ESTIMATE Analytic conjugate posterior for the stored sample.
            arguments
                obj (1,1) minnesotamniwbvarm
                Y   = obj.Y
                opts.Display
                opts.X
                opts.Y0
            end

            obj.assertStoredSample(Y, "estimate");
            args = namedargs2cell(opts);
            [MN, Summary] = estimate@conjugatebvarm(obj, obj.Y, args{:});

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
            %SIMULATE Draw coefficients and covariance given the stored sample.
            arguments
                obj (1,1) minnesotamniwbvarm
                Y   = obj.Y
                opts.NumDraws
                opts.X
                opts.Y0
            end
            obj.assertStoredSample(Y, "simulate");
            args = namedargs2cell(opts);
            [varargout{1:nargout}] = simulate@conjugatebvarm(obj, obj.Y, args{:});
        end

        function varargout = forecast(obj, numperiods, Y, varargin)
            %FORECAST Forecast responses beyond the stored sample.
            arguments
                obj        (1,1) minnesotamniwbvarm
                numperiods (1,1) double {mustBeInteger, mustBePositive}
                Y          = obj.Y
            end
            arguments (Repeating)
                varargin
            end
            obj.assertStoredSample(Y, "forecast");
            [varargout{1:nargout}] = ...
                forecast@conjugatebvarm(obj, numperiods, obj.Y, varargin{:});
        end

        function varargout = simsmooth(obj, Y, varargin)
            %SIMSMOOTH Simulation smoother, defaulting to the stored sample.
            %   Deliberately NOT sample-checked. simsmooth is the primitive
            %   the inherited machinery calls internally - bvar/forecast in
            %   particular passes the stored sample extended with NaN rows
            %   over the forecast horizon - so a strict equality check here
            %   would reject the toolbox's own legitimate calls. The check
            %   belongs on the user-facing entry points above.
            arguments
                obj (1,1) minnesotamniwbvarm
                Y   = obj.Y
            end
            arguments (Repeating)
                varargin
            end
            [varargout{1:nargout}] = ...
                simsmooth@conjugatebvarm(obj, Y, varargin{:});
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
            [Mu, V, Omega, DoF] = minnesotamniwbvarm.updateNIW( ...
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

        function psi = resolvePsi(obj, Psi, Y, numlags)
            %RESOLVEPSI Residual-variance scale from a vector or an estimator name.
            if isstring(Psi) || ischar(Psi)
                method = string(Psi);
                mustBeMember(method, ["exact","conditional"]);
                psi = estimateResidualVariances(Y, numlags, Method=method);
                return
            end

            if ~isnumeric(Psi) || ~isvector(Psi) || numel(Psi) ~= obj.NumSeries
                error("minnesotamniwbvarm:invalidPsi", ...
                    "Psi must be a 1-by-%d numeric vector or one of " + ...
                    """exact"" / ""conditional"".", obj.NumSeries);
            end
            if any(~isfinite(Psi)) || any(Psi <= 0)
                error("minnesotamniwbvarm:invalidPsi", ...
                    "Psi values must be finite and positive.");
            end
            psi = reshape(double(Psi), 1, []);
        end

        function assertStoredSample(obj, Y, caller)
            %ASSERTSTOREDSAMPLE Reject a sample other than the one built from.
            if istabular(Y)
                Y = Y{:,:};
            end
            if ~isequal(Y, obj.Y)
                error("minnesotamniwbvarm:sampleMismatch", ...
                    "%s was called with a sample that differs from the one " + ...
                    "this prior was built from. Psi and the lambda4/lambda5 " + ...
                    "dummy observations derive from the stored sample, so a " + ...
                    "different sample needs a new minnesotamniwbvarm.", caller);
            end
        end

    end

    % ---- shared conjugate NIW algebra -----------------------------------
    methods (Static, Access = private)

        function [Y, nvp2] = resolveSample(Y, numseries, numlags, nvp2)
            %RESOLVESAMPLE Validate Y and adopt tabular variable names.
            if istabular(Y)
                if ~isfield(nvp2, "SeriesNames")
                    nvp2.SeriesNames = string(Y.Properties.VariableNames);
                end
                Y = Y{:,:};
            end

            if ~isnumeric(Y) || ~ismatrix(Y)
                error("minnesotamniwbvarm:invalidData", ...
                    "Y must be a numeric matrix or a table/timetable.");
            end
            if size(Y, 2) ~= numseries
                error("minnesotamniwbvarm:invalidData", ...
                    "Y must have %d columns, one per series.", numseries);
            end
            if size(Y, 1) <= numlags
                error("minnesotamniwbvarm:invalidData", ...
                    "Y must have more rows than the lag order P = %d.", numlags);
            end
            Y = double(Y);
        end

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
