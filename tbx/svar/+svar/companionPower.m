function Phi = companionPower(mdl, maxLag)
%COMPANIONPOWER  Selected numSeries-by-numSeries blocks of companion
%matrix powers 0..maxLag.
%
%   Phi = COMPANIONPOWER(mdl, maxLag) returns an
%   numSeries-by-numSeries-by-(maxLag+1) array with
%   Phi(:,:,h+1) = selector*comp^h*selector', h = 0,...,maxLag, where comp
%   is companionMatrix(mdl) and selector picks out the top numSeries rows
%   and columns.
%
%   This is the expensive part shared by svar.irf, svar.fevd, and
%   svar.identification.ShockPath: it depends only on mdl, never on the
%   impact matrix or shock column, so callers that need it for more than
%   one impact/shock (e.g. svar.Identification.identify, looping over
%   shocks within a draw) should compute it once per mdl and slice/contract
%   the result themselves rather than recomputing it per shock.

arguments
    mdl
    maxLag (1, 1) double {mustBeInteger}
end

numSeries = mdl.NumSeries;
comp = svar.companionMatrix(mdl);
stateDim = size(comp, 1);
selector = [eye(numSeries) zeros(numSeries, stateDim - numSeries)];

Phi = zeros(numSeries, numSeries, maxLag + 1);
Apow = eye(stateDim);
for h = 0:maxLag
    Phi(:, :, h + 1) = selector * Apow * selector';
    Apow = comp * Apow;
end
end
