function spec = minnesotaSpec(Method, varargin)

arguments
    Method (1,1) string = "mniw"
end

arguments (Repeating)
    varargin
end

validMethods = ["mniw","inw","normal"];
method = lower(Method);
specArgs = varargin;

if ~any(method == validMethods)
    if isempty(varargin)
        error("minnesotaSpec:unknownMethod", ...
            "Method must be one of: %s.", strjoin(validMethods, ", "));
    end
    method = "mniw";
    specArgs = [{Method}, varargin];
end

switch method
    case "mniw"
        spec = minnesotamniwSpec(specArgs{:});
    case "inw"
        spec = minnesotainwSpec(specArgs{:});
    case "normal"
        spec = minnesotanSpec(specArgs{:});
end

end
