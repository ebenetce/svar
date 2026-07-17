function spec = minnesotaSpec(varargin)

[aliasMap, aliasSummary] = localAliasMap();

if nargin == 0
    Method = "mniw";
    specArgs = {};
else
    firstArg = varargin{1};
    if (isstring(firstArg) || ischar(firstArg)) && isscalar(string(firstArg))
        firstMethod = lower(erase(string(firstArg), ["-","_"," "]));
        if isKey(aliasMap, firstMethod) || nargin == 1
            Method = string(firstArg);
            specArgs = varargin(2:end);
        else
            Method = "mniw";
            specArgs = varargin;
        end
    else
        Method = "mniw";
        specArgs = varargin;
    end
end

method = lower(erase(string(Method), ["-","_"," "]));
if ~isKey(aliasMap, method)
    error("minnesotaSpec:unknownMethod", ...
        "Unrecognized Method '%s'.%s%s", Method, newline, aliasSummary);
end
method = aliasMap(method);

switch method
    case "mniw"
        spec = minnesotamniwSpec.create(specArgs{:});
    case "inw"
        spec = minnesotainwSpec.create(specArgs{:});
    case "normal"
        spec = minnesotanSpec.create(specArgs{:});
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
