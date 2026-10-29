classdef SignRestriction < svar.identification.Restriction
    %SIGNRESTRICTION - Restrict an impulse response to one sign
    %   R = SIGNRESTRICTION(SIGN,Response=RESPONSE,Shock=SHOCK,Horizon=HORIZON)
    %   requires the named RESPONSE to the named SHOCK to have the
    %   specified sign at every integer HORIZON.
    %
    %   SIGNRESTRICTION properties:
    %       Response - Name of the restricted response variable
    %       Shock    - Name of the restricted structural shock
    %       Horizon  - Nonnegative impulse-response horizons
    %       Sign     - Required sign, either 1 or -1
    %
    %   See also svar.identification.RestrictionSet, svar.identify

    properties (SetAccess = immutable)
        Response (1,1) string  % Restricted response-variable name
        Horizon (1,:) double   % Nonnegative impulse-response horizons
        Sign (1,1) double      % Required response sign
    end

    properties (SetAccess = private, Hidden)
        ResponseIndex (1,1) double = NaN
        ShockIndex (1,1) double = NaN
    end

    methods
        function obj = SignRestriction(sign, nvp)
            %SIGNRESTRICTION - Create an impulse-response sign restriction
            arguments
                sign (1,1) double {mustBeMember(sign, [-1 1])}
                nvp.Response (1,1) string
                nvp.Shock (1,1) string
                nvp.Horizon (1,:) double {mustBeInteger, mustBeNonnegative}
            end

            if strlength(nvp.Response) == 0 || strlength(nvp.Shock) == 0
                error("svar:SignRestriction:EmptyName", ...
                    "Response and Shock must be nonempty names.");
            end

            obj.Response = nvp.Response;
            obj.Shock = nvp.Shock;
            obj.Horizon = unique(nvp.Horizon, "stable");
            obj.Sign = sign;
        end

        function [value, lower, upper] = evaluate(obj, state)
            %EVALUATE - Evaluate signed impulse responses for a draw
            value = obj.Sign*reshape( ...
                state.IRF(obj.Horizon + 1, obj.ResponseIndex, obj.ShockIndex), ...
                [], 1);
            lower = zeros(size(value));
            upper = inf(size(value));
        end

        function obj = resolve(obj, seriesNames, shockNames)
            %RESOLVE - Cache response and shock indices for evaluation
            [hasResponse, obj.ResponseIndex] = ismember(obj.Response, seriesNames);
            if ~hasResponse
                error("svar:SignRestriction:UnknownResponse", ...
                    "Response ""%s"" is not defined for this identification.", ...
                    obj.Response);
            end

            [hasShock, obj.ShockIndex] = ismember(obj.Shock, shockNames);
            if ~hasShock
                error("svar:SignRestriction:UnknownShock", ...
                    "Shock ""%s"" is not defined for this identification.", ...
                    obj.Shock);
            end
        end

        function horizon = requiredHorizon(obj)
            %REQUIREDHORIZON - Return the largest restricted horizon
            horizon = max(obj.Horizon);
        end
    end
end
