classdef hyperprior
    %HYPERPRIOR A single hyperprior on one scalar hyperparameter.
    %
    %   A lightweight value object holding ONE distribution family and its
    %   canonical coefficients, plus a stable log-density and the pieces an
    %   optimiser needs (box bounds, a starting point). It is the third state
    %   of a Minnesota spec field, alongside the existing scalar-or-bounds
    %   convention:
    %
    %       spec.lambda1 = 0.2;                          % FIXED
    %       spec.lambda1 = [0.05 0.5];                   % FREE, flat prior
    %       spec.lambda1 = hyperprior("Gamma", 2, 0.1);  % FREE, Gamma prior
    %
    %   so the prior on a hyperparameter IS the value of that hyperparameter -
    %   there is no separate name-keyed hyperprior container to keep in sync,
    %   and a misspelt parameter name fails as a bad property assignment
    %   rather than silently dropping a prior term.
    %
    %   Parameterisations
    %   -----------------
    %   The constructor takes NATIVE parameters (always closed-form, no
    %   solver, safe on an optimiser hot path):
    %       hyperprior("Gamma",        k,     theta)   % shape, scale
    %       hyperprior("InverseGamma", alpha, beta)    % shape, scale
    %       hyperprior("Beta",         alpha, beta)
    %
    %   The FROMMOMENTS factory takes an interpretable (mode, sd) pair and
    %   converts ONCE, here, at definition:
    %       hyperprior.fromMoments("Gamma", mode, sd)  % closed form
    %       hyperprior.fromMoments("Beta",  mode, sd)  % 1-D fzero (base MATLAB)
    %       hyperprior.fromMoments("InverseGamma", mode, sd)  % 1-D fzero
    %   No Optimization Toolbox, no eqnproblem: Gamma is analytic, Beta and
    %   Inverse-Gamma each reduce to a single scalar root solved with FZERO.
    %
    %   Methods
    %   -------
    %       logpdf(hp, x)   stable log density (log space throughout; returns
    %                       -Inf outside the support rather than throwing, so
    %                       an optimiser can probe freely)
    %       bounds(hp)      [lo hi] quantile box for the optimiser
    %       x0(hp)          starting point (the mode, or the median as fallback)
    %       modeOf/meanOf   convenience moments
    %
    %   Dependencies: GAMINV/BETAINV (Statistics and Machine Learning Toolbox)
    %   are used by BOUNDS only, and are already required transitively by the
    %   Econometrics Toolbox this extension targets. LOGPDF depends on base
    %   MATLAB only (gammaln, betaln).
    %
    %   See also minnesotaBaseSpec, minnesotaSpec.

    properties (SetAccess = private)
        Distribution (1,1) string
        Params       (1,2) double    % [shape scale] (Gamma/IG) or [alpha beta] (Beta)
        Bounds       (1,2) double     
        X0           (1,1) double
    end

    methods

        function obj = hyperprior(distribution, mode, sd, nvp)
            %HYPERPRIOR Construct from native parameters (no solver).
            arguments
                distribution (1,1) string {mustBeMember(distribution, ["Gamma", "InverseGamma", "Beta"])}
                mode (1,1) double {mustBePositive}
                sd (1,1) double {mustBePositive}
                nvp.Bounds (1,2) double
                nvp.X0 (1,1) double
            end  

            obj.Distribution = distribution;
            [p1, p2] = fromMoments(distribution, mode, sd);
            obj.Params       = [p1, p2];
                
            if ~isfield(nvp, "Bounds")                
                obj.Bounds = obj.quantileBounds;
            else
                if nvp.Bounds(1) >= nvp.Bounds(2)
                    error("hyperpriors:BoundsWrongOrder", "Lower Bounds can't be equal or larger than the Upper Bound")
                end
                obj.Bounds = nvp.Bounds;
            end

            if ~isfield(nvp, "X0")
                obj.X0 = obj.initialValue;
            else
                obj.X0 = nvp.X0;
            end
        end

        function lp = logpdf(obj, x)
            %LOGPDF Log density at x. Vectorised; -Inf outside the support.
            arguments
                obj (1,1) hyperprior
                x   double
            end
            a = obj.Params(1);
            b = obj.Params(2);
            lp = -inf(size(x));

            switch obj.Distribution
                case "Gamma"                              % [k theta]
                    ok = x > 0;
                    lp(ok) = (a - 1).*log(x(ok)) - x(ok)./b ...
                             - a.*log(b) - gammaln(a);

                case "InverseGamma"                       % [alpha beta]
                    ok = x > 0;
                    lp(ok) = a.*log(b) - gammaln(a) ...
                             - (a + 1).*log(x(ok)) - b./x(ok);

                case "Beta"                               % [alpha beta]
                    ok = x > 0 & x < 1;
                    lp(ok) = (a - 1).*log(x(ok)) + (b - 1).*log(1 - x(ok)) ...
                             - betaln(a, b);
            end
        end

        function [lo, hi] = quantileBounds(obj, massBounds)
            %BOUNDS Quantile box [lo hi] for an optimiser.
            %   Default mass bounds [0.001 0.999]. For a free hyperparameter
            %   carrying this prior, these are the box constraints FMINCON
            %   searches within.
            arguments
                obj        (1,1) hyperprior
                massBounds (1,2) double {mustBeBetween(massBounds, 0, 1)} = [1e-3, 1-1e-3]
            end
            a = obj.Params(1);
            b = obj.Params(2);
            switch obj.Distribution
                case "Gamma"
                    lo = gaminv(massBounds(1), a, b);
                    hi = gaminv(massBounds(2), a, b);
                case "InverseGamma"
                    % X ~ IG(a,b)  <=>  1/X ~ Gamma(a, 1/b).
                    lo = 1 / gaminv(massBounds(2), a, 1/b);
                    hi = 1 / gaminv(massBounds(1), a, 1/b);
                case "Beta"
                    lo = betainv(massBounds(1), a, b);
                    hi = betainv(massBounds(2), a, b);
            end
            if nargout <= 1
                lo = [lo hi];   % allow b = hp.bounds()
            end
        end

        function x = initialValue(obj)
            %X0 Optimiser starting point: the mode where defined, else median.
            m = obj.modeOf();
            if isfinite(m) && m > 0 && (obj.Distribution ~= "Beta" || m < 1)
                x = m;
                return
            end
            a = obj.Params(1);
            b = obj.Params(2);
            switch obj.Distribution                       % median via icdf(0.5)
                case "Gamma",        x = gaminv(0.5, a, b);
                case "InverseGamma", x = 1 / gaminv(0.5, a, 1/b);
                case "Beta",         x = betainv(0.5, a, b);
            end
        end

        function m = modeOf(obj)
            %MODEOF Distribution mode (may be a boundary value / undefined).
            a = obj.Params(1);
            b = obj.Params(2);
            switch obj.Distribution
                case "Gamma"                              % (k-1)*theta, k>1
                    if a > 1, m = (a - 1)*b; else, m = 0; end
                case "InverseGamma"
                    m = b / (a + 1);
                case "Beta"                               % (a-1)/(a+b-2), a,b>1
                    if a > 1 && b > 1
                        m = (a - 1)/(a + b - 2);
                    else
                        m = NaN;
                    end
            end
        end

        function mu = meanOf(obj)
            %MEANOF Distribution mean (Inf/undefined where it does not exist).
            a = obj.Params(1);
            b = obj.Params(2);
            switch obj.Distribution
                case "Gamma",        mu = a*b;
                case "InverseGamma"                       % beta/(alpha-1), a>1
                    if a > 1, mu = b/(a - 1); else, mu = Inf; end
                case "Beta",         mu = a/(a + b);
            end
        end

    end

end

function [p1, p2] = fromMoments(distribution, mode, sd)
%FROMMOMENTS Build from an interpretable (mode, sd) pair.
%   Gamma is analytic; Beta and Inverse-Gamma each reduce to one
%   scalar root (FZERO). Conversion happens HERE, once, so the
%   returned object holds native params and never re-solves.
arguments
    distribution (1,1) string ...
        {mustBeMember(distribution, ["Gamma","InverseGamma","Beta"])}
    mode (1,1) double {mustBePositive}
    sd   (1,1) double {mustBePositive}
end
switch distribution
    case "Gamma"
        [p1, p2] = gammaFromMoments(mode, sd);
    case "InverseGamma"
        [p1, p2] = invGammaFromMoments(mode, sd);
    case "Beta"
        [p1, p2] = betaFromMoments(mode, sd);
end
end

function [k, theta] = gammaFromMoments(m, s)
%GAMMAFROMMOMENTS Closed form. mode=(k-1)*theta, var=k*theta^2.
%   With u = sqrt(k): m/s = u - 1/u  =>  u^2 - (m/s)u - 1 = 0.
r     = m/s;
u     = (r + sqrt(r^2 + 4))/2;      % positive root; u>1 for r>0
k     = u^2;                        % => k>1, interior mode
theta = s/u;                        % = s/sqrt(k)
end

function [alpha, beta] = invGammaFromMoments(m, s)
%INVGAMMAFROMMOMENTS mode=beta/(alpha+1); one root alpha>2 always.
%   (s/m)^2 = (alpha+1)^2 / ((alpha-1)^2 (alpha-2)).
target = (s/m)^2;
f = @(al) (al + 1).^2 ./ ((al - 1).^2 .* (al - 2)) - target;
alpha = solveScalar(f, 2 + 1e-9, 1e6, ...
    "InverseGamma", m, s);
beta  = m*(alpha + 1);
end

function [alpha, beta] = betaFromMoments(m, s)
%BETAFROMMOMENTS Reduce to concentration c=alpha+beta, root c>2.
%   a=m(c-2)+1, b=c-a, var=a*b/(c^2 (c+1)). Feasible sd is bounded
%   above at a given mode (max at c=2, Beta(1,1)); an infeasible
%   (mode, sd) is reported rather than silently mis-solved.
g = @(c) betaVarAtConcentration(c, m) - s^2;
c = solveScalar(g, 2 + 1e-9, 1e6, "Beta", m, s);
alpha = m*(c - 2) + 1;
beta  = c - alpha;
end

function v = betaVarAtConcentration(c, m)
a = m.*(c - 2) + 1;
b = c - a;
v = a.*b ./ (c.^2 .* (c + 1));
end

function root = solveScalar(f, lo, hi, name, m, s)
%SOLVESCALAR Bracketed 1-D root with a clear feasibility error.
flo = f(lo);
fhi = f(hi);
if ~isfinite(flo) || ~isfinite(fhi) || sign(flo) == sign(fhi)
    error("hyperprior:infeasibleMoments", ...
        "No %s matches mode=%g, sd=%g (the requested spread is " + ...
        "outside the feasible range for this mode).", name, m, s);
end
root = fzero(f, [lo hi]);
end