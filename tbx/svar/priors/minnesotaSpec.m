function spec = minnesotaSpec(Method, varargin)

arguments
    Method (1,1) string
end

arguments (Repeating)
    varargin
end

[aliasMap, aliasSummary] = localAliasMap();

method = lower(erase(string(Method), ["-","_"," "]));
if ~isKey(aliasMap, method)
    error("minnesotaSpec:unknownMethod", ...
        "Unrecognized Method '%s'.%s%s", Method, newline, aliasSummary);
end
method = aliasMap(method);

switch method
    case "mniw"
        spec = minnesotamniwSpec.create(varargin{:});
    case "inw"
        spec = minnesotainwSpec.create(varargin{:});
    case "normal"
        spec = minnesotanSpec.create(varargin{:});
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