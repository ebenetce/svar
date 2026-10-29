classdef (Abstract) Restriction < matlab.mixin.Heterogeneous
    %RESTRICTION - Abstract structural-identification restriction
    %   Restrictions evaluate a structural-draw state and define a closed
    %   interval of admissible values. svar.identify combines public restriction
    %   declarations into an internal RestrictionSet.
    %
    %   See also svar.SignRestriction, svar.identification.RestrictionSet

    properties
        Shock (1,1) string  % Named structural shock
    end

    methods (Abstract)
        [value, lower, upper] = evaluate(obj, state)
    end

    methods
        function tf = isSatisfied(obj, state)
            %ISSATISFIED - Test whether a state satisfies this restriction
            [value, lower, upper] = obj.evaluate(state);
            tf = all(isfinite(value), "all") ...
                && all(value >= lower, "all") ...
                && all(value <= upper, "all");
        end

        function horizon = requiredHorizon(~)
            %REQUIREDHORIZON - Return the largest required IRF horizon
            horizon = 0;
        end
    end
end
