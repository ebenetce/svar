function [irf, Phi] = irf(varMdl, impact, horizon, nvp)

arguments
    varMdl 
    impact 
    horizon 
    nvp.Phi double = [];
end

% assert size(Phi) == [numSeries, numSeries, horizon + 1],

numSeries = size(impact, 1);
numShocks = size(impact, 2);

if isempty(nvp.Phi)
    Phi = svar.companionPower(varMdl, horizon);
else
    Phi = nvp.Phi;
end

irf = zeros(horizon + 1, numSeries, numShocks);
for h = 0:horizon
    irf(h + 1, :, :) = Phi(:, :, h + 1) * impact;
end

end