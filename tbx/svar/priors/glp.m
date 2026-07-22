function [mdl, info] = glp(numseries, numlags, Y, Psi, nvp, nvp2)
%GLP Tune Minnesota hyperparameters by marginal likelihood.
%   mdl = GLP(numseries,numlags,Y,Psi,Spec=spec) tunes the free fields in
%   SPEC and, optionally, Psi by maximising the analytic marginal likelihood
%   of a MINNESOTAMNIWBVARM prior against Y.
%
%   SPEC fields are fixed when scalar, free with a flat prior when given as
%   [lower upper], and free with an embedded log prior when given as a
%   HYPERPRIOR object.
%
%   Psi is one of:
%       1-by-n numeric          fixed residual variances
%       2-by-n numeric          free residual variances with flat priors
%       scalar hyperprior       broadcast independently to all n variances
%       1-by-n hyperprior       free residual variances with log priors
%
%   Hyperprior objects provide Bounds, X0, and logpdf. There is no separate
%   PriorCoef or PsiBand input.
%
%   Example
%   -------
%       spec = minnesotaSpec("mniw", ...
%           lambda1=hyperprior("Gamma",0.2,0.4,Bounds=[1e-4 5]), ...
%           lambda4=hyperprior("Gamma",1,1,Bounds=[1e-4 50]), ...
%           lambda5=hyperprior("Gamma",1,1,Bounds=[1e-4 50]));
%       psi0 = estimateResidualVariances(Y,1,Method="conditional");
%       Psi = arrayfun(@(x) hyperprior("InverseGamma",0.02^2,0.02^2, ...
%           X0=x, Bounds=[1/100 100]*x), psi0);
%       mdl = glp(size(Y,2),4,Y,Psi,Spec=spec);
%
%   See also MINNESOTASPEC, MINNESOTAMNIWBVARM, HYPERPRIOR.

arguments
    numseries (1,1) double {mustBeInteger, mustBePositive}
    numlags   (1,1) double {mustBeInteger, mustBePositive}
    Y         {mustBeNonempty}
    Psi
    nvp.Spec          (1,1) svar.minnesotaBaseSpec = minnesotaSpec("mniw")
    nvp.OptimOptions      = optimoptions("fmincon", ...
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

if ~isa(nvp.Spec, "svar.minnesotamniwSpec")
    error("glp:unsupportedSpec", ...
        "GLP requires a minnesotaSpec(""mniw"") because it uses an analytic marginal likelihood.");
end

if ~istabular(Y) && ~isnumeric(Y)
    error("glp:BadDataType", "Responses input must be numeric or tabular")
end

if istabular(Y)
    if ~isfield(nvp2, "SeriesNames")
        nvp2.SeriesNames = string(Y.Properties.VariableNames);
    end
    Y = Y{:,:};
end

if size(Y, 2) ~= numseries
    error("glp:invalidData", ...
        "Y must have %d columns, one per series.", numseries);
end

buildArgs = namedargs2cell(nvp2);

baseSpec = nvp.Spec;
[lamX0, lamLB, lamUB, lambdaNames] = baseSpec.pack();
[psiX0, psiLB, psiUB, psiFixed, psiPrior, psiNames] = localPackPsi(Psi, numseries);

x0 = [lamX0, psiX0];
lb = [lamLB, psiLB];
ub = [lamUB, psiUB];
nLambdaFree = numel(lambdaNames);
psiFree = ~isempty(psiX0);
useHyperprior = ~isempty(baseSpec.hyperpriorFields()) || ~isempty(psiPrior);

    function [candSpec, candPsi] = unpackAll(x)
        candSpec = baseSpec.unpack(x(1:nLambdaFree), lambdaNames);
        if psiFree
            candPsi = x(nLambdaFree+1:end);
        else
            candPsi = psiFixed;
        end
    end

    function negObj = objective(x)
        [candSpec, candPsi] = unpackAll(x);
        try
            mdlCandidate = candSpec.build(numseries, numlags, candPsi, buildArgs{:});
            negObj = -logMarginalLikelihood(mdlCandidate, Y);
        catch
            negObj = 1e10;
        end

        if isfinite(negObj) && negObj <= 1e10
            logPrior = baseSpec.logHyperprior(x(1:nLambdaFree), lambdaNames) ...
                + localLogPsiHyperprior(psiPrior, candPsi);
            negObj = negObj - logPrior;
        end

        if ~isfinite(negObj) || negObj > 1e10
            negObj = 1e10;
        end
    end

if isempty(x0)
    xHat = x0;
    fminconOutput = struct("Skipped", true, ...
        "Reason", "No free hyperparameters.");
    fval = objective(xHat);
    exitflag = NaN;
else
    [xHat, fval, exitflag, fminconOutput] = fmincon(@objective, x0, ...
        [], [], [], [], lb, ub, [], nvp.OptimOptions);
end

[finalSpec, finalPsi] = unpackAll(xHat);
mdl = finalSpec.build(numseries, numlags, finalPsi, buildArgs{:});

if nargout > 1
    info = struct( ...
        "LambdaNames",     lambdaNames, ...
        "PsiNames",        psiNames, ...
        "PsiFree",         psiFree, ...
        "InitialSpec",     baseSpec, ...
        "FinalSpec",       finalSpec, ...
        "InitialPsi",      Psi, ...
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

function [x0, lb, ub, fixedPsi, psiPrior, names] = localPackPsi(Psi, numseries)
names = "Psi" + string(1:numseries);
fixedPsi = [];
psiPrior = [];

if isa(Psi, "hyperprior")
    psiPrior = reshape(Psi, 1, []);
    if isscalar(psiPrior) && numseries > 1
        psiPrior = repmat(psiPrior, 1, numseries);
    end
    if numel(psiPrior) ~= numseries
        error("glp:invalidPsi", ...
            "Psi hyperprior input must be scalar or have one element per series.");
    end
    x0 = arrayfun(@(p) p.X0, psiPrior);
    lb = arrayfun(@(p) p.Bounds(1), psiPrior);
    ub = arrayfun(@(p) p.Bounds(2), psiPrior);
    return
end

if ~isnumeric(Psi)
    error("glp:invalidPsi", ...
        "Psi must be numeric or a hyperprior array.");
end

if isvector(Psi) && numel(Psi) == numseries
    fixedPsi = reshape(Psi, 1, []);
    if any(~isfinite(fixedPsi)) || any(fixedPsi <= 0)
        error("glp:invalidPsi", ...
            "Fixed Psi values must be finite and positive.");
    end
    x0 = [];
    lb = [];
    ub = [];
    return
end

if ismatrix(Psi) && isequal(size(Psi), [2 numseries])
    lb = Psi(1,:);
    ub = Psi(2,:);
    if any(~isfinite(lb)) || any(~isfinite(ub)) || any(lb <= 0) || any(lb >= ub)
        error("glp:invalidPsi", ...
            "Psi bounds must be finite, positive, and ordered lower < upper.");
    end
    x0 = sqrt(lb .* ub);
    return
end

error("glp:invalidPsi", ...
    "Psi must be 1-by-n numeric, 2-by-n numeric bounds, scalar hyperprior, or 1-by-n hyperprior.");
end

function lp = localLogPsiHyperprior(psiPrior, psi)
lp = 0;
if isempty(psiPrior)
    return
end
for i = 1:numel(psi)
    lp = lp + psiPrior(i).logpdf(psi(i));
end
end
