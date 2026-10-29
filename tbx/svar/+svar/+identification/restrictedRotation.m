function rotation = restrictedRotation(baseImpact, restrictions)
%RESTRICTEDROTATION - Draw a rotation conditional on impact zero restrictions

arguments
    baseImpact (:,:) double {mustBeReal, mustBeFinite}
    restrictions (1,1) svar.identification.RestrictionSet
end

numSeries = size(baseImpact, 1);
if ~restrictions.hasZeroRestrictions()
    rotation = svar.utils.haarRotation(numSeries);
    return
end

rotation = zeros(numSeries);
zeroRestricted = false(1, numSeries);
for shock = 1:numSeries
    zeroRestricted(shock) = ~isempty(restrictions.zeroRows(shock, baseImpact));
end
drawOrder = [find(zeroRestricted), find(~zeroRestricted)];
placedColumns = zeros(numSeries);
numPlaced = 0;

for shock = drawOrder
    zeroRows = restrictions.zeroRows(shock, baseImpact);
    nullBasis = null([zeroRows; placedColumns(:, 1:numPlaced).']);
    if isempty(nullBasis)
        error("svar:identify:OverIdentifiedZeroRestrictions", ...
            "Shock ""%s"" has no remaining direction after its zero " + ...
            "restrictions and previously drawn shocks.", ...
            restrictions.ShockOrder(shock));
    end

    direction = randn(size(nullBasis, 2), 1);
    rotation(:, shock) = nullBasis*(direction/norm(direction));
    numPlaced = numPlaced + 1;
    placedColumns(:, numPlaced) = rotation(:, shock);
end
end
