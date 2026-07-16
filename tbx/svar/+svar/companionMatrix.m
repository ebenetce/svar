function companion = companionMatrix(mdl)

arguments
    mdl
end

numSeries = mdl.NumSeries;
numLags = mdl.P;

companion = zeros(numSeries * numLags);
companion(1:numSeries, :) = cat(2, mdl.AR{:});

if numLags > 1
    companion(numSeries + 1:end, 1:end - numSeries) = ...
        eye(numSeries * (numLags - 1));
end

end