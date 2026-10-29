classdef RestrictionSet
    %RESTRICTIONSET - Combine structural-identification restrictions
    %   SET = RESTRICTIONSET(SERIESNAMES) creates an internal restriction
    %   accumulator for SERIESNAMES. svar.identify constructs this object
    %   from public restriction declarations.
    %
    %   RESTRICTIONSET functions:
    %       isSatisfied    - Test a structural draw against every restriction
    %       requiredHorizon - Return the largest required IRF horizon
    %
    %   RESTRICTIONSET properties:
    %       Restrictions - Cell array of constituent restriction objects
    %
    %   See also svar.SignRestriction, svar.ZeroRestriction, svar.identify

    properties (SetAccess = private)
        SeriesNames (1,:) string  % Model response-variable names
        ShockOrder (1,:) string   % Complete structural-shock order
        Restrictions (1,:) cell   % Constituent restriction objects
    end

    methods
        function obj = RestrictionSet(seriesNames, nvp)
            %RESTRICTIONSET - Create an internal restriction accumulator
            arguments
                seriesNames (1,:) string
                nvp.ShockOrder (1,:) string
            end

            if numel(unique(seriesNames)) ~= numel(seriesNames) || ...
                    numel(unique(nvp.ShockOrder)) ~= numel(nvp.ShockOrder)
                error("svar:RestrictionSet:DuplicateNames", ...
                    "SeriesNames and ShockOrder must contain unique names.");
            end

            obj.SeriesNames = seriesNames;
            obj.ShockOrder = nvp.ShockOrder;
            obj.Restrictions = cell(1, 0);
        end

        function obj = add(obj, restriction)
            %ADD - Add a public restriction declaration
            if ~isa(restriction, "svar.identification.Restriction")
                error("svar:RestrictionSet:InvalidRestriction", ...
                    "Restriction must inherit from svar.identification.Restriction.");
            end
            if ~any(obj.ShockOrder == restriction.Shock)
                error("svar:RestrictionSet:UnknownShock", ...
                    "Shock ""%s"" is not in ShockOrder.", restriction.Shock);
            end
            if isa(restriction, "svar.SignRestriction") || ...
                    isa(restriction, "svar.ZeroRestriction")
                restriction = restriction.resolve( ...
                    obj.SeriesNames, obj.ShockOrder);
            end
            obj.Restrictions{end + 1} = restriction;
        end

        function tf = hasZeroRestrictions(obj)
            %HASZERORESTRICTIONS - Return true when impact zeros are present
            tf = any(cellfun(@(r) isa(r, "svar.ZeroRestriction"), ...
                obj.Restrictions));
        end

        function rows = zeroRows(obj, shockIndex, baseImpact)
            %ZEROROWS - Return rotation rows that a shock must null
            arguments
                obj
                shockIndex (1,1) double {mustBeInteger, mustBePositive}
                baseImpact (:,:) double {mustBeReal, mustBeFinite}
            end

            zeroRestrictions = obj.Restrictions(cellfun( ...
                @(r) isa(r, "svar.ZeroRestriction") && ...
                r.appliesToShock(shockIndex), obj.Restrictions));
            rows = zeros(numel(zeroRestrictions), size(baseImpact, 2));
            for restrictionIndex = 1:numel(zeroRestrictions)
                rows(restrictionIndex, :) = ...
                    zeroRestrictions{restrictionIndex}.rotationRow(baseImpact);
            end
        end

        function selector = zeroSelectorMatrix(obj, shockIndex)
            %ZEROSELECTORMATRIX - Select impact responses fixed at zero
            arguments
                obj
                shockIndex (1,1) double {mustBeInteger, mustBePositive}
            end

            zeroRestrictions = obj.Restrictions(cellfun( ...
                @(r) isa(r, "svar.ZeroRestriction") && ...
                r.appliesToShock(shockIndex), obj.Restrictions));
            selector = zeros(numel(zeroRestrictions), numel(obj.SeriesNames));
            for restrictionIndex = 1:numel(zeroRestrictions)
                selector(restrictionIndex, :) = ...
                    zeroRestrictions{restrictionIndex}.responseSelector( ...
                    numel(obj.SeriesNames));
            end
        end

        function tf = isSatisfied(obj, state)
            %ISSATISFIED - Test a structural draw against every restriction
            tf = true;
            for restriction = obj.Restrictions
                if ~restriction{1}.isSatisfied(state)
                    tf = false;
                    return
                end
            end
        end

        function values = margins(obj, state)
            %MARGINS - Return finite signed distances from restriction bounds
            values = zeros(0, 1);
            for restriction = obj.Restrictions
                [value, lower, upper] = restriction{1}.evaluate(state);
                value = value(:);
                lower = lower(:);
                upper = upper(:);
                values = [values; value(isfinite(lower)) - ...
                    lower(isfinite(lower)); upper(isfinite(upper)) - ...
                    value(isfinite(upper))]; %#ok<AGROW>
            end
        end

        function horizon = requiredHorizon(obj)
            %REQUIREDHORIZON - Return the largest required IRF horizon
            horizons = cellfun(@requiredHorizon, obj.Restrictions);
            horizon = max(horizons);
        end
    end
end
