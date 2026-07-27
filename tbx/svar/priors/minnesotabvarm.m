function mdl = minnesotabvarm(numseries, numlags, Y, varargin)
%MINNESOTABVARM Minnesota prior for a Bayesian VAR.
%   PriorMdl = MINNESOTABVARM(NUMSERIES,NUMLAGS,Y) creates a Minnesota
%   (Lit
% terman) prior for the sample Y. Everything data-dependent - the
%   residual-variance scale, and any dummy observations - is resolved here,
%   so the object that comes back is a complete prior.
%
%   PriorMdl = MINNESOTABVARM(...,Type=TYPE) selects the prior family:
%
%     TYPE                 Class                     Posterior      Extra
%     -------------------- ------------------------- -------------- ----------
%     "mniw" (default)     svar.minnesotamniwbvarm   analytic NIW   lambda4/5
%     "inw"                svar.minnesotainwbvarm    Gibbs          lambda2
%     "normal"             svar.minnesotanbvarm      analytic       lambda2
%
%   Aliases, ignoring hyphens, underscores, and spaces:
%     mniw    conjugate, matrixnormal
%     inw     independent, kadiyala, kadiyalakarlsson, semiconjugate
%     normal  litterman, fixed, fixedsigma
%
%   PriorMdl = MINNESOTABVARM(...,Psi=PSI) sets the residual-variance scale:
%   a 1-by-NUMSERIES vector, or "exact" (default) / "conditional" naming the
%   estimator ESTIMATERESIDUALVARIANCES applies to Y. Pass PSI numerically
%   inside a tuning loop - the string form refits an AR per series.
%
%   PriorMdl = MINNESOTABVARM(...,Name=Value) sets the hyperparameters of the
%   chosen family (lambda1, lambda3, Vc, PriorMean for all; lambda4/lambda5
%   for "mniw"; lambda2 for "inw" and "normal") and any option the underlying
%   BVAR accepts (SeriesNames, IncludeConstant, IncludeTrend, NumPredictors,
%   Description). Options are validated by the family's own constructor, so
%   passing lambda2 to "mniw" - where conjugacy pins it to 1 - is an error,
%   not a silent no-op.
%
%   Example
%   -------
%       prior = minnesotabvarm(3, 4, Y);                        % conjugate
%       prior = minnesotabvarm(3, 4, Y, lambda4=1, lambda5=1);  % + dummies
%       prior = minnesotabvarm(3, 4, Y, Type="inw", lambda2=0.5);
%       post  = estimate(prior);                                % sample is stored
%
%   NOTE. This is a dispatching function, not a class: a MATLAB constructor
%   must return an object of its own class, so a single classdef cannot both
%   be the shared base of the three families and hand back one of them. The
%   common ancestor is svar.minnesotabvarmBase - use ISA against that to test
%   for "some Minnesota prior". This mirrors BAYESVARM in Econometrics
%   Toolbox, which likewise dispatches on a model-type argument.
%
%   See also SVAR.MINNESOTAMNIWBVARM, SVAR.MINNESOTAINWBVARM,
%   SVAR.MINNESOTANBVARM, GLP, ESTIMATERESIDUALVARIANCES.

arguments
    numseries (1,1) double {mustBeInteger, mustBePositive}
    numlags   (1,1) double {mustBeInteger, mustBePositive}
    Y                {mustBeNonempty}
end

arguments (Repeating)
    varargin
end

[type, args] = localExtractType(varargin);

switch type
    case "mniw"
        mdl = svar.minnesotamniwbvarm(numseries, numlags, Y, args{:});
    case "inw"
        mdl = svar.minnesotainwbvarm(numseries, numlags, Y, args{:});
    case "normal"
        mdl = svar.minnesotanbvarm(numseries, numlags, Y, args{:});
end

end

function [type, args] = localExtractType(args)
%LOCALEXTRACTTYPE Pull Type out of the name-value list and resolve its alias.
type = "mniw";

names = args(1:2:end);
isType = cellfun(@(n) (isstring(n) || ischar(n)) && strcmpi(n, "Type"), names);

if ~any(isType)
    return
end

if nnz(isType) > 1
    error("minnesotabvarm:repeatedType", "Type was specified more than once.");
end

position = 2*find(isType) - 1;
if position + 1 > numel(args)
    error("minnesotabvarm:missingTypeValue", "Type requires a value.");
end

type = localResolveAlias(args{position + 1});
args([position, position + 1]) = [];
end

function type = localResolveAlias(value)
if ~(isstring(value) || ischar(value)) || ~isscalar(string(value))
    error("minnesotabvarm:invalidType", "Type must be a string scalar.");
end

[aliasMap, aliasSummary] = localAliasMap();
key = lower(erase(string(value), ["-","_"," "]));

if ~isKey(aliasMap, key)
    error("minnesotabvarm:unknownType", ...
        "Unrecognized Type '%s'.%s%s", value, newline, aliasSummary);
end

type = aliasMap(key);
end

function [aliasMap, aliasSummary] = localAliasMap()
mniwAliases   = ["mniw","conjugate","matrixnormal"];
inwAliases    = ["inw","independent","kadiyala","kadiyalakarlsson","semiconjugate"];
normalAliases = ["normal","litterman","fixed","fixedsigma"];

aliasMap = dictionary(mniwAliases, "mniw", inwAliases, "inw", ...
    normalAliases, "normal");

aliasSummary = strjoin([ ...
    "Available Type aliases:", ...
    "  mniw   (Matrix-Normal-Inverse-Wishart): " + strjoin(mniwAliases, ", "), ...
    "  inw    (Independent Normal-Wishart):    " + strjoin(inwAliases, ", "), ...
    "  normal (fixed-Sigma Normal):            " + strjoin(normalAliases, ", "), ...
    "Hyphens, underscores, and spaces are ignored when matching aliases."], ...
    newline);
end
