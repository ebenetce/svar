function [mdl, info] = glp(numseries, numlags, Y, Psi, nvp, nvp2)
%GLP Tune Minnesota hyperparameters by (log) marginal likelihood.
%   mdl = GLP(numseries, numlags, Y, Psi) tunes whichever fields
%   of NVP.SPEC are left FREE (2-element [lower upper] bounds) by maximising
%   the analytic marginal likelihood of a MINNESOTABVARM prior against Y, via
%   FMINCON, and returns the tuned MINNESOTABVARM. Fields left as a scalar in
%   NVP.SPEC are held fixed. See MINNESOTASPEC for the scalar-or-bounds
%   convention.
%
%   [mdl, info] = GLP(...) also returns a struct with the
%   optimisation details (bounds, starting point, FMINCON output, psi used).
%
%   Name-value arguments
%   ---------------------
%   Spec          A MINNESOTASPEC. Any field left as a scalar is fixed; any
%                 field given as [lower upper] is optimised. Default
%                 minnesotaSpec() (everything fixed at its default -> no
%                 lambda is tuned unless you pass a Spec with some field
%                 given as a range).
%   PsiBand       Scalar (default 1) -> psi is FIXED at PsiBand * baseline,
%                 where baseline is the per-series AR(PsiLags) residual
%                 variance (see PsiLags/PsiMethod below). A 2-element
%                 [lower upper] -> psi is FREE, one variable per series, each
%                 bounded to [lower, upper] * that series' baseline - the
%                 same multiplicative-band convention used historically for
%                 this search (e.g. [1/100, 100]).
%   PriorCoef     Struct passed to MINNESOTASPEC.logHyperprior at each
%                 candidate point. Default empty struct -> pure marginal-
%                 likelihood (ML) tuning. Supply fields to switch to MAP
%                 (posterior-mode) tuning - see MINNESOTASPEC.logHyperprior.
%   IncludeConstant, IncludeTrend, SeriesNames
%                 Forwarded to MINNESOTASPEC.build / MINNESOTABVARM.
%                 NumPredictors is not supported (the marginal-likelihood /
%                 dummy path does not carry exogenous regressors).
%   OptimOptions  Options for FMINCON. Default a tight-tolerance interior-
%                 point-family configuration suitable for this smooth,
%                 analytic objective.
%
%   Any objective evaluation that fails or is non-finite (e.g. a candidate
%   point makes the dummy-augmented prior numerically degenerate) is clipped
%   to a large finite penalty (1e10) rather than propagated as Inf/NaN, so
%   FMINCON's line search is never handed a non-finite value.
%
%   psi is computed from Y ONCE up front (either fixed via PsiBand, or given
%   bounds to search) and held or searched independently of the lambdas -
%   matching the standard GLP algorithm, where psi is a data-derived scale,
%   not something whose baseline is refit per candidate. Each objective
%   evaluation is a cheap MINNESOTASPEC.build (no refitting) plus one
%   analytic marginal likelihood - no sampler, no rng needed.
%
%   Example
%   -------
%       spec = minnesotaSpec(lambda1 = [1e-3, 5], lambda4 = [1e-3, 5], ...
%                             lambda5 = [1e-3, 5]);
%       [mdl, info] = GLP(3, 4, Y, Psi, ...
%           Spec = spec, PsiBand = [1/100, 100], IncludeConstant = true);
%       Posterior = mdl.estimate(Y);
%
%   See also MINNESOTASPEC, MINNESOTABVARM, ESTIMATERESIDUALVARIANCES.

arguments
    numseries (1,1) double {mustBeInteger, mustBePositive}
    numlags   (1,1) double {mustBeInteger, mustBePositive}
    Y         {mustBeNonempty}
    Psi
    nvp.Spec          (1,1) minnesotamniwSpec = minnesotamniwSpec()
    nvp.PsiBand       (1,:) double {mustBeScalarOrPositiveBounds} = 1    
    nvp.PriorCoef     (1,1) struct = struct()
    nvp.OptimOptions    = optimoptions("fmincon", ...
        Display = "final", ...
        FiniteDifferenceStepSize = 1e-4, ...
        FunctionTolerance = 1e-12, ...
        StepTolerance = 1e-12, ...
        ConstraintTolerance = 1e-12)
    nvp2.Description
    nvp2.IncludeConstant
    nvp2.IncludeTrend
    nvp2.NumPredictors
    nvp2.SeriesNames
end

if ~istabular(Y) && ~isnumeric(Y)
    error('glp:BadDataType', 'Responses input must be numeric or tabular')
end

if istabular(Y)
    if ~isfield(nvp, "SeriesNames")
        nvp.SeriesNames = string(Y.Properties.VariableNames);
    end
    Y = Y{:,:};     
end

if size(Y, 2) ~= numseries
    error("glp:invalidData", ...
        "Y must have %d columns, one per series.", numseries);
end

% ---- build-time options ---------------------------------------------
buildArgs = namedargs2cell(nvp2);

useHyperprior = ~isempty(fieldnames(nvp.PriorCoef));

% ---- assemble the combined [free lambdas ; free psi] decision vector -----
baseSpec         = nvp.Spec;
[lamX0, lamLB, lamUB, lambdaNames] = baseSpec.pack();

psiFree = numel(nvp.PsiBand) == 2;
if psiFree
    psiLB = Psi * nvp.PsiBand(1);
    psiUB = Psi * nvp.PsiBand(2);
    psiX0 = sqrt(psiLB .* psiUB);          % geometric-mean start, per series
else
    psiLB = []; psiUB = []; psiX0 = [];
end

x0 = [lamX0, psiX0];
lb = [lamLB, psiLB];
ub = [lamUB, psiUB];
nLambdaFree = numel(lambdaNames);

    function [candSpec, candPsi] = unpackAll(x)
        candSpec = baseSpec.unpack(x(1:nLambdaFree), lambdaNames);
        if psiFree
            candPsi = x(nLambdaFree+1:end);
        else
            candPsi = Psi * nvp.PsiBand;   % PsiBand is a fixed scalar here
        end
    end

    function negObj = objective(x)
        [candSpec, candPsi] = unpackAll(x);
        try
            mdlCandidate = candSpec.build(numseries, numlags, candPsi, buildArgs{:});
            negObj = mdlCandidate.negativeLogMarginalLikelihood(Y);
        catch
            negObj = 1e10;
        end
        if ~isfinite(negObj) || negObj > 1e10
            negObj = 1e10;
        end
        if useHyperprior
            negObj = negObj - candSpec.logHyperprior(nvp.PriorCoef, candPsi);
        end
    end

% ---- optimise (skip FMINCON entirely if nothing is free) -----------------
if isempty(x0)
    xHat = x0;
    fminconOutput = struct("Skipped", true, ...
        "Reason", "No free hyperparameters (Spec fully fixed and PsiBand scalar).");
    fval = objective(xHat); % kept for info even though there was nothing to search
    exitflag = NaN;
else
    [xHat, fval, exitflag, fminconOutput] = fmincon(@objective, x0, ...
        [], [], [], [], lb, ub, [], nvp.OptimOptions);
end

[finalSpec, finalPsi] = unpackAll(xHat);
mdl = finalSpec.build(numseries, numlags, finalPsi, buildArgs{:});

if nargin > 1
info = struct( ...
    "LambdaNames",     lambdaNames, ...
    "PsiFree",         psiFree, ...
    "InitialSpec",     baseSpec, ...
    "FinalSpec",       finalSpec, ...
    "BaselinePsi",     Psi, ...
    "FinalPsi",        finalPsi, ...
    "UsedHyperprior",  useHyperprior, ...
    "X0",              x0, ...
    "LowerBound",      lb, ...
    "UpperBound",      ub, ...
    "XHat",            xHat, ...
    "Objective",       fval, ...
    "ExitFlag",        exitflag, ...
    "FminconOutput",   fminconOutput);
end

end

function mustBeScalarOrPositiveBounds(x)
%MUSTBESCALARORPOSITIVEBOUNDS Validator for PsiBand: a positive scalar
%   (fixed multiplier) or a finite 2-element [lower upper] with lower < upper
%   (free multiplicative band).
if numel(x) ~= 1 && numel(x) ~= 2
    error("glp:invalidPsiBand", ...
        "PsiBand must be a scalar or a 2-element [lower upper]; got %d elements.", ...
        numel(x));
end
if any(x <= 0) || any(~isfinite(x))
    error("glp:invalidPsiBand", ...
        "PsiBand must be finite and positive.");
end
if numel(x) == 2 && x(1) >= x(2)
    error("glp:invalidPsiBand", ...
        "PsiBand bounds must satisfy lower < upper.");
end
end
