classdef ZeroRestriction < svar.identification.Restriction
    %ZERORESTRICTION - Restrict an impact response to zero
    %   R = ZERORESTRICTION(Response=RESPONSE,Shock=SHOCK) requires the
    %   named RESPONSE to have zero impact response to the named SHOCK.
    %
    %   Zero restrictions are constructed by the identification sampler,
    %   rather than tested by rejection sampling.
    %
    %   See also svar.SignRestriction, svar.identify

    properties (SetAccess = immutable)
        Response (1,1) string  % Restricted response-variable name
        Horizon (1,1) double   % Restricted impulse-response horizon
    end

    properties (SetAccess = private, Hidden)
        ResponseIndex (1,1) double = NaN
        ShockIndex (1,1) double = NaN
    end

    methods
        function obj = ZeroRestriction(nvp)
            %ZERORESTRICTION - Create an impact zero restriction
            arguments
                nvp.Response (1,1) string
                nvp.Shock (1,1) string
                nvp.Horizon (1,1) double {mustBeInteger, mustBeNonnegative} = 0
            end

            if strlength(nvp.Response) == 0 || strlength(nvp.Shock) == 0
                error("svar:ZeroRestriction:EmptyName", ...
                    "Response and Shock must be nonempty names.");
            end
            if nvp.Horizon ~= 0
                error("svar:ZeroRestriction:UnsupportedHorizon", ...
                    "ZeroRestriction currently supports only Horizon=0.");
            end

            obj.Response = nvp.Response;
            obj.Shock = nvp.Shock;
            obj.Horizon = nvp.Horizon;
        end

        function [value, lower, upper] = evaluate(obj, state)
            %EVALUATE - Evaluate the restricted impact response
            value = state.Impact(obj.ResponseIndex, obj.ShockIndex);
            tolerance = 100*eps(max(1, norm(state.Impact, "fro")));
            lower = -tolerance;
            upper = tolerance;
        end

        function obj = resolve(obj, seriesNames, shockNames)
            %RESOLVE - Cache response and shock indices for evaluation
            [hasResponse, obj.ResponseIndex] = ismember(obj.Response, seriesNames);
            if ~hasResponse
                error("svar:ZeroRestriction:UnknownResponse", ...
                    "Response ""%s"" is not defined for this identification.", ...
                    obj.Response);
            end

            [hasShock, obj.ShockIndex] = ismember(obj.Shock, shockNames);
            if ~hasShock
                error("svar:ZeroRestriction:UnknownShock", ...
                    "Shock ""%s"" is not defined for this identification.", ...
                    obj.Shock);
            end
        end

        function tf = appliesToShock(obj, shockIndex)
            %APPLIESTOSHOCK - Test whether this restriction targets a shock
            tf = obj.ShockIndex == shockIndex;
        end

        function row = rotationRow(obj, baseImpact)
            %ROTATIONROW - Return the rotation row constrained to zero
            row = baseImpact(obj.ResponseIndex, :);
        end

        function selector = responseSelector(obj, numSeries)
            %RESPONSESELECTOR - Select the restricted response from impact
            selector = zeros(1, numSeries);
            selector(obj.ResponseIndex) = 1;
        end
    end
end
