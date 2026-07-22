function mdl = varmFromCoefficients(template, coefficients, covariance)
%VARMFROMCOEFFICIENTS Reconstruct a fitted VARM from a BVAR draw.
%   mdl = VARMFROMCOEFFICIENTS(template, coefficients, covariance) builds a
%   concrete VARM model from one draw of AR/constant/trend/predictor
%   coefficients and one draw of the innovations covariance, using TEMPLATE
%   (any CONJUGATEBVARM - or a subclass such as MINNESOTABVARM, prior or
%   posterior) only for its structural layout: NumSeries, P, SeriesNames,
%   IncludeConstant, IncludeTrend, NumPredictors.
%
%   COEFFICIENTS is an m-by-NumSeries matrix in the same row layout as
%   TEMPLATE.Mu reshaped to m-by-NumSeries - i.e. lag blocks first
%   (P*NumSeries rows), then the constant row (if IncludeConstant), then the
%   trend row (if IncludeTrend), then predictor rows (if NumPredictors > 0).
%   This is exactly the shape returned per-draw by SIMULATE on a
%   CONJUGATEBVARM (or MINNESOTABVARM).
%
%   This function is prior-agnostic: it has no Minnesota-specific logic and
%   belongs to neither MINNESOTABVARM nor MINNESOTASPEC. Typical use is
%   turning simulate() draws into VARM objects for per-draw FORECAST / IRF /
%   FEVD (credible bands), independent of which prior produced the draw.
%
%   See also CONJUGATEBVARm, SIMULATE, VARM.

arguments
    template bvar
    coefficients (:,:) double {mustBeReal, mustBeFinite}
    covariance   (:,:) double {mustBeReal, mustBeFinite}
end

numSeries = template.NumSeries;
numLags   = template.P;

if template.NumPredictors > 0
    error("varmFromCoefficients:predictorsUnsupported", ...
        "VARM has no exogenous-predictor slot analogous to Beta here; " + ...
        "this reconstruction does not support NumPredictors > 0.");
end

expectedRows = numLags*numSeries + template.IncludeConstant + template.IncludeTrend;
if size(coefficients, 1) ~= expectedRows || size(coefficients, 2) ~= numSeries
    error("varmFromCoefficients:sizeMismatch", ...
        "coefficients must be %d-by-%d for this template (got %d-by-%d).", ...
        expectedRows, numSeries, size(coefficients, 1), size(coefficients, 2));
end

mdl             = varm(numSeries, numLags);
mdl.SeriesNames = template.SeriesNames;

for lag = 1:numLags
    rows        = ((lag - 1)*numSeries + 1):(lag*numSeries);
    mdl.AR{lag} = coefficients(rows, :).';
end

row = numLags*numSeries;

if template.IncludeConstant
    row = row + 1;
    mdl.Constant = coefficients(row, :).';
else
    mdl.Constant = zeros(numSeries, 1);
end

if template.IncludeTrend
    row = row + 1;
    mdl.Trend = coefficients(row, :).';
end

mdl.Covariance = covariance;

end