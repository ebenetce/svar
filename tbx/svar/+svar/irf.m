function [irf, Phi] = irf(varMdl, impact, horizon, nvp)
%IRF - Compute impulse responses from a VAR model and impact matrix
%   RESPONSES = IRF(varMdl,IMPACT,HORIZON) computes impulse responses
%   from horizon 0 through HORIZON for each shock column in IMPACT.
%
%   RESPONSES = IRF(...,Phi=PHI) uses precomputed moving-average
%   coefficient blocks. Use this form when reusing PHI across apply calls
%   for the same VAR model and horizon.
%
%   [RESPONSES,PHI] = IRF(...) also returns the moving-average
%   coefficient blocks used in the calculation.
%
%   See also fevd, companionPower, varm

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
