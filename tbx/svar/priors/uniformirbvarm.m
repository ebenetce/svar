classdef uniformirbvarm < conjugatebvarm
    %uniformirbvarm Reduced-form prior induced by uniform IR parameters.
    %
    % This class implements the reduced-form component of the prior in
    % Arias, Rubio-Ramirez, and Waggoner, "Uniform Priors for Impulse
    % Responses." It excludes the structural orthogonal-matrix block. For a
    % VAR with n variables and k regressors per equation, the reduced-form
    % prior density is proportional to det(Sigma)^((k + s)/2), where the
    % paper's impulse-response prior uses s = -3.
    %
    % The density is represented as the equivalent improper conjugate NIW
    % limit: diffuse coefficients, zero inverse-Wishart scale, and
    % DoF = -2*k - n - 1 - s. MATLAB's inherited estimate method then
    % produces the paper's reduced-form NIW posterior.

    properties (SetAccess = private)
        DeterminantShift (1,1) double {mustBeReal, mustBeFinite} = -3
    end

    properties (Dependent, Hidden, SetAccess = private)
        NumEquationCoefficients
        LogDetExponent
    end

    methods
        function obj = uniformirbvarm(numseries, numlags, nvp, nvp2)
            arguments
                numseries 
                numlags 
                nvp.?conjugatebvarm
                nvp2.DeterminantShift
            end

            args = namedargs2cell(nvp);
            obj@conjugatebvarm(numseries, numlags, args{:});

            if isfield(nvp2, "DeterminantShift")
                obj.DeterminantShift = nvp2.DeterminantShift;
            end

            k = obj.NumEquationCoefficients;

            obj.Mu = zeros(numseries*k, 1);
            obj.V = diag(inf(k, 1));
            obj.Omega = zeros(numseries, numseries);

            warningState = warning("query", "econ:bvar:bvar:TooSmallDoF");
            cleanup = onCleanup(@() warning(warningState.state, ...
                "econ:bvar:bvar:TooSmallDoF"));
            warning("off", "econ:bvar:bvar:TooSmallDoF");

            obj.DoF = obj.priorDoF();
        end

        function PostMdl = estimate(obj, varargin)
            % Estimate model
            PostMdlbase = estimate@conjugatebvarm(obj, varargin{:});

            % Create posterior
            PostMdl = conjugatebvarm(obj.NumSeries, obj.P, ...
                IncludeConstant = obj.IncludeConstant, ...
                IncludeTrend = obj.IncludeTrend, ...
                NumPredictors = obj.NumPredictors, ...
                SeriesNames = PostMdlbase.SeriesNames);
            
            PostMdl.PrivateMu = PostMdlbase.Mu;
            PostMdl.PrivateV = PostMdlbase.V;
            PostMdl.PrivateOmega = PostMdlbase.Omega;
            PostMdl.PrivateDoF = PostMdlbase.DoF;
        end

        function value = get.NumEquationCoefficients(obj)
            value = obj.NumSeries*obj.P ...
                + double(obj.IncludeConstant) ...
                + double(obj.IncludeTrend) ...
                + obj.NumPredictors;
        end

        function value = get.LogDetExponent(obj)
            value = obj.logDetExponent();
        end
       
    end

    methods (Access = private)

        function dof = priorDoF(obj)
            dof = -2*obj.NumEquationCoefficients ...
                - obj.NumSeries - 1 - obj.DeterminantShift;
        end

        function exponent = logDetExponent(obj)
            exponent = (obj.NumEquationCoefficients + obj.DeterminantShift)/2;
        end

    end

end
