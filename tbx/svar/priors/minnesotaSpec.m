function spec = minnesotaSpec(Method, varargin)

arguments
    Method (1,1) string = "mniw"
end

arguments (Repeating)
    varargin
end

switch Method
    case "mniw"
        spec = minnesotamniwSpec(varargin{:});
    case "inw"
        spec = minnesotaiwSpec(varargin{:});
    case "normal"
        spec = minnesotanSpec(varargin{:});
end

end