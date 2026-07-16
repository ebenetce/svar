function psi = estimateResidualVariances(Y, numLags, opts)
%ESTIMATERESIDUALVARIANCES Per-series residual variance from a univariate AR(P) fit.
%   psi = ESTIMATERESIDUALVARIANCES(Y, numLags) fits an independent AR(numLags)
%   to each column of Y and returns each fit's innovations variance as a
%   1-by-NumSeries row vector - the standard "ppsi" scale consumed by
%   MINNESOTASPEC.build / MINNESOTABVARM.
%
%   psi = ESTIMATERESIDUALVARIANCES(..., Method="conditional") uses conditional
%   maximum likelihood - closed-form, equivalent to OLS on the observed sample,
%   currently implemented via VARM/ESTIMATE - instead of the default "exact"
%   method (approximately exact Gaussian ML, currently implemented via
%   ARIMA/ESTIMATE, which handles the first P observations through the
%   process's stationary/backcast distribution rather than discarding them).
%   The two are NOT numerically identical, especially in short samples -
%   "exact" is the default because it matches the estimator used historically
%   in this codebase; "conditional" is available because it is a direct,
%   non-iterative solve and noticeably faster when NumSeries is large and
%   exact reproducibility with legacy results is not required.
%
%   Method values name the STATISTICAL choice, not the MATLAB function that
%   implements it - e.g. "conditional" is currently VARM-based, but that
%   implementation detail can change without changing this argument's meaning
%   or any caller's code.
%
%   This is a plain, prior-agnostic utility - it belongs to neither
%   MINNESOTABVARM nor MINNESOTASPEC. Compute it ONCE and reuse the result
%   across repeated builds (e.g. inside a hyperparameter search) rather than
%   recomputing it per candidate.
%
%   See also MINNESOTASPEC, MINNESOTABVARM, ARIMA, VARM.

arguments
    Y
    numLags (1,1) double {mustBeInteger, mustBePositive}
    opts.Method    (1,1) string {mustBeMember(opts.Method, ["exact","conditional"])} = "exact"
end

numSeries = size(Y, 2);

if istabular(Y)
    Y = Y{:,:};
end

psi       = zeros(1, numSeries);

for series = 1:numSeries
    y = Y(:, series);
    switch opts.Method
        case "exact"
            Mdl         = arima(numLags, 0, 0);
            EstMdl      = estimate(Mdl, y, Display="off");
            psi(series) = EstMdl.Variance;
        case "conditional"
            Mdl         = varm(1, numLags);
            EstMdl      = estimate(Mdl, y, Display="off");
            psi(series) = EstMdl.Covariance;
    end
end

end
