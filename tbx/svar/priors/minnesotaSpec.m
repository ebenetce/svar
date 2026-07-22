function spec = minnesotaSpec(prior, varargin)

arguments
    prior (1,1) string
end

arguments (Repeating)
    varargin
end

[aliasMap, aliasSummary] = localAliasMap();

prior = lower(erase(string(prior), ["-","_"," "]));
if ~isKey(aliasMap, prior)
    error("minnesotaSpec:unknownMethod", ...
        "Unrecognized Method '%s'.%s%s", prior, newline, aliasSummary);
end
prior = aliasMap(prior);

switch prior
    case "mniw"
        spec = svar.minnesotamniwSpec.create(varargin{:});
    case "inw"
        spec = svar.minnesotainwSpec.create(varargin{:});
    case "normal"
        spec = svar.minnesotanSpec.create(varargin{:});
end

end

function [aliasMap, aliasSummary] = localAliasMap()
mniwAliases = ["mniw","conjugate","matrixnormal"];
inwAliases = ["inw","independent","kadiyala","kadiyalakarlsson","semiconjugate"];
normalAliases = ["normal","litterman","fixed","fixedsigma"];

aliasMap = dictionary(mniwAliases, "mniw", inwAliases, "inw", normalAliases, "normal");

aliasSummary = strjoin([ ...
    "Available Method aliases:", ...
    "  mniw   (Matrix-Normal-Inverse-Wishart): " + strjoin(mniwAliases, ", "), ...
    "  inw    (Independent Normal-Wishart):    " + strjoin(inwAliases, ", "), ...
    "  normal (fixed-Sigma Normal):            " + strjoin(normalAliases, ", "), ...
    "Hyphens, underscores, and spaces are ignored when matching aliases."], ...
    newline);

end
