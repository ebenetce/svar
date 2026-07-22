function logML = logMarginalLikelihood(Mdl, Y, opts)
%LOGMARGINALLIKELIHOOD Log marginal likelihood for supported BVAR priors.
%   LOGML = LOGMARGINALLIKELIHOOD(MDL,Y) returns log p(Y | MDL) for proper
%   conjugate BVAR and fixed-Sigma Normal BVAR prior objects.
%
%   LOGML = LOGMARGINALLIKELIHOOD(MDL) uses the sample stored on MDL, for
%   prior objects that carry one (see MINNESOTAMNIWBVARM).
%
%   LOGML = LOGMARGINALLIKELIHOOD(MDL,Y,X=X,Y0=Y0) supplies exogenous
%   predictors or presample responses to the sufficient-statistics
%   calculation.
%
%   [LOGML,DETAILS] = LOGMARGINALLIKELIHOOD(___) also returns posterior
%   moments used by the analytic calculation.
%
%   Semiconjugate and improper priors are not supported by this function.
%
%   See also MARGINALLIKELIHOOD, CONJUGATEBVARM, NORMALBVARM

arguments
    Mdl (1,1)
    Y = []
    opts.X = []
    opts.Y0 = []
end

Y = localResolveSample(Mdl, Y);
Y = localNumericData(Y, "Y");
X = localNumericData(opts.X, "X");
Y0 = localNumericData(opts.Y0, "Y0");

if isa(Mdl, "semiconjugatebvarm")
    error("logMarginalLikelihood:unsupportedModel", ...
        "Analytic log marginal likelihood is not implemented for %s.", class(Mdl));
end

if isa(Mdl, "diffusebvarm") || isa(Mdl, "weakbvarm") ...
        || isa(Mdl, "uniformirbvarm")
    error("logMarginalLikelihood:improperPrior", ...
        "Log marginal likelihood is not defined for the improper prior %s.", ...
        class(Mdl));
end

if isa(Mdl, "normalbvarm")
    [XX, XY, YY, numObs] = localSufficientStatistics(Mdl, Y, X, Y0);
    logML = localNormalLogML(Mdl, XX, XY, YY, numObs);
    return
end

% NB: no Minnesota special case. MINNESOTAMNIWBVARM materialises its
% lambda4/lambda5 dummy observations in the constructor, so by the time a
% prior reaches here its Mu/V/Omega/DoF are already the augmented ones and
% the generic conjugate path is correct.
if isa(Mdl, "conjugatebvarm")
    [XX, XY, YY, numObs] = localSufficientStatistics(Mdl, Y, X, Y0);
    logML = localConjugateLogML(Mdl, XX, XY, YY, numObs);
    return
end

error("logMarginalLikelihood:unsupportedModel", ...
    "Log marginal likelihood is not implemented for %s.", class(Mdl));
end

function Y = localResolveSample(Mdl, Y)
%LOCALRESOLVESAMPLE Fall back to a sample carried by the model object.
if ~isempty(Y)
    return
end

if isprop(Mdl, "Y") && ~isempty(Mdl.Y)
    Y = Mdl.Y;
    return
end

error("logMarginalLikelihood:missingData", ...
    "Y is required because %s does not carry a sample.", class(Mdl));
end

function A = localNumericData(A, name)
if isempty(A)
    return
end

if istabular(A)
    A = A{:,:};
end

if ~isnumeric(A)
    error("logMarginalLikelihood:invalidData", ...
        "%s must be numeric or tabular.", name);
end
end

function [XX, XY, YY, numObs] = localSufficientStatistics(Mdl, Y, X, Y0)
if size(Y, 2) ~= Mdl.NumSeries
    error("logMarginalLikelihood:invalidData", ...
        "Y must have %d columns, one per series.", Mdl.NumSeries);
end

if Mdl.NumPredictors > 0 && isempty(X)
    error("logMarginalLikelihood:missingPredictors", ...
        "X is required because the model has NumPredictors = %d.", ...
        Mdl.NumPredictors);
end

if Mdl.NumPredictors == 0 && ~isempty(X)
    error("logMarginalLikelihood:unexpectedPredictors", ...
        "X was supplied but the model has NumPredictors = 0.");
end

[XX, XY, YY, numObs] = bvar.sufficientStatistics(Y, Y0, X, ...
    Mdl.NumSeries, Mdl.P, Mdl.IncludeConstant, Mdl.IncludeTrend, ...
    Mdl.NumPredictors, false);
end

function logML = localConjugateLogML(Mdl, XX, XY, YY, numObs)
n = Mdl.NumSeries;
V = Mdl.V;
Omega = Mdl.Omega;
DoF = Mdl.DoF;

localAssertProperConjugate(V, Omega, DoF, n, class(Mdl));

k = size(V, 1);
priorDoF = DoF;
postDoF = DoF + numObs;

try
    B0 = reshape(Mdl.Mu, k, n);
    V0inv = V \ eye(k);
    prec = V0inv + XX;
    Bn = prec \ (V0inv*B0 + XY);

    resid = YY - Bn'*XY - XY'*Bn + Bn'*XX*Bn;
    resid = (resid + resid')/2;
    shift = Bn - B0;
    incr = resid + shift'*(V0inv*shift);
    incr = (incr + incr')/2;

    Lo = chol((Omega + Omega')/2, "lower");
    bbb = Lo \ incr / Lo';
    eb = real(eig((bbb + bbb')/2));
    eb(eb < 0) = 0;
    sumOmegaRatio = sum(log(eb + 1));
    logDetOmega0 = 2*sum(log(diag(Lo)));

    D = chol((V + V')/2, "lower");
    aaa = D'*XX*D;
    ea = real(eig((aaa + aaa')/2));
    ea(ea < 0) = 0;
    sumVRatio = sum(log(ea + 1));

    logML = -0.5*numObs*n*log(pi) ...
        + localLogMvGamma(0.5*postDoF, n) ...
        - localLogMvGamma(0.5*priorDoF, n) ...
        - 0.5*numObs*logDetOmega0 ...
        - 0.5*postDoF*sumOmegaRatio ...
        - 0.5*n*sumVRatio;
    
catch
    logML = -Inf;    
end
end

function logML = localNormalLogML(Mdl, XX, XY, YY, numObs)
n = Mdl.NumSeries;
Sigma = (Mdl.Sigma + Mdl.Sigma')/2;
priorCovariance = (Mdl.V + Mdl.V')/2;

localAssertProperNormal(priorCovariance, Sigma, class(Mdl));

try
    sigmaInv = Sigma \ eye(n);
    priorPrecision = priorCovariance \ eye(size(priorCovariance, 1));
    posteriorPrecision = priorPrecision + kron(sigmaInv, XX);
    posteriorV = posteriorPrecision \ eye(size(posteriorPrecision, 1));
    posteriorV = (posteriorV + posteriorV')/2;

    dataMoment = XY*sigmaInv;
    posteriorMu = posteriorV*(priorPrecision*Mdl.Mu(:) + dataMoment(:));

    sigmaFactor = chol(Sigma, "lower");
    priorFactor = chol(priorCovariance, "lower");
    posteriorPrecisionFactor = chol( ...
        (posteriorPrecision + posteriorPrecision')/2, "lower");

    logDetSigma = 2*sum(log(diag(sigmaFactor)));
    logDetPrior = 2*sum(log(diag(priorFactor)));
    logDetPosteriorPrecision = ...
        2*sum(log(diag(posteriorPrecisionFactor)));

    priorMean = Mdl.Mu(:);
    dataQuadratic = trace(sigmaInv*YY);
    priorQuadratic = priorMean'*(priorPrecision*priorMean);
    posteriorQuadratic = posteriorMu'*(posteriorPrecision*posteriorMu);

    logML = -0.5*(numObs*n*log(2*pi) ...
        + numObs*logDetSigma ...
        + logDetPrior ...
        + logDetPosteriorPrecision ...
        + dataQuadratic ...
        + priorQuadratic ...
        - posteriorQuadratic);
   
catch
    logML = -Inf;
end
end

function localAssertProperConjugate(V, Omega, DoF, n, modelClass)
if any(~isfinite(V), "all") || any(~isfinite(Omega), "all") ...
        || ~isfinite(DoF) || DoF <= n - 1
    error("logMarginalLikelihood:improperPrior", ...
        "Log marginal likelihood requires a proper conjugate prior; %s is improper.", ...
        modelClass);
end
end

function localAssertProperNormal(V, Sigma, modelClass)
if any(~isfinite(V), "all") || any(~isfinite(Sigma), "all")
    error("logMarginalLikelihood:improperPrior", ...
        "Log marginal likelihood requires a proper Normal prior; %s is improper.", ...
        modelClass);
end
end

function value = localLogMvGamma(a, dimension)
j = 1:dimension;
value = dimension*(dimension - 1)*0.25*log(pi) ...
    + sum(gammaln(a + 0.5*(1 - j)));
end