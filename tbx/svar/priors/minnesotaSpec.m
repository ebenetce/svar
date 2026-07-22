function spec = minnesotaSpec(prior, varargin)
%MINNESOTASPEC Hyperparameter recipe for a Minnesota BVAR prior.
%
%   TRANSITIONAL. The conjugate ("mniw") family no longer goes through a
%   spec: construct it directly with MINNESOTAMNIWBVARM(numseries,numlags,Y,
%   ...), and tune it with GLP, which now takes the hyperparameters itself.
%   The "inw" and "normal" families still use this factory until they are
%   converted to the same constructor-takes-data API.
%
%   See also MINNESOTAMNIWBVARM, GLP.

arguments
    prior (1,1) string
end

arguments (Repeating)
    varargin
end

[aliasMap, aliasSummary] = localAliasMap();

prior = lower(erase(string(prior), ["-","_"," "]));

if isKey(localRetiredAliasMap(), prior)
    error("minnesotaSpec:retiredMethod", ...
        "The conjugate Minnesota prior no longer uses a spec. " + ...
        "Call minnesotamniwbvarm(numseries, numlags, Y, ...) directly, " + ...
        "or glp(numseries, numlags, Y, ...) to tune it.");
end

if ~isKey(aliasMap, prior)
    error("minnesotaSpec:unknownMethod", ...
        "Unrecognized Method '%s'.%s%s", prior, newline, aliasSummary);
end
prior = aliasMap(prior);

switch prior
    case "inw"
        spec = svar.minnesotainwSpec.create(varargin{:});
    case "normal"
        spec = svar.minnesotanSpec.create(varargin{:});
end

end

function map = localRetiredAliasMap()
map = dictionary(["mniw","conjugate","matrixnormal"], true);
end

function [aliasMap, aliasSummary] = localAliasMap()
inwAliases = ["inw","independent","kadiyala","kadiyalakarlsson","semiconjugate"];
normalAliases = ["normal","litterman","fixed","fixedsigma"];

aliasMap = dictionary(inwAliases, "inw", normalAliases, "normal");

aliasSummary = strjoin([ ...
    "Available Method aliases:", ...
    "  inw    (Independent Normal-Wishart):    " + strjoin(inwAliases, ", "), ...
    "  normal (fixed-Sigma Normal):            " + strjoin(normalAliases, ", "), ...
    "The conjugate (mniw) prior is constructed directly - see MINNESOTAMNIWBVARM.", ...
    "Hyphens, underscores, and spaces are ignored when matching aliases."], ...
    newline);

end
