function [decomposition, responses, Phi] = fevd(varMdl, impact, horizon, nvp)
%FEVD - Forecast error variance decomposition from structural responses
%   decomposition = FEVD(varMdl, impact, horizon) computes the share of each
%   variable's forecast error variance attributable to each shock column in
%   impact, from horizon 0 through horizon.
%
%   decomposition = FEVD(..., Phi=Phi) uses precomputed moving-average
%   coefficient blocks. Use this form when reusing Phi across apply calls for
%   the same VAR model and horizon.
%
%   [decomposition, responses, Phi] = FEVD(...) also returns the impulse
%   responses and moving-average coefficient blocks used in the calculation.
%
%   See also irf, companionPower, varm

arguments
    varMdl
    impact (:,:) double {mustBeReal, mustBeFinite}
    horizon (1,1) double {mustBeInteger, mustBeNonnegative}
    nvp.Phi double = []
end

args = namedargs2cell(nvp);
[responses, Phi] = svar.irf(varMdl, impact, horizon, args{:});

cumulativeSquaredResponses = cumsum(responses.^2, 1);
forecastErrorVariance = sum(cumulativeSquaredResponses, 3);
decomposition = cumulativeSquaredResponses ./ forecastErrorVariance;

end
