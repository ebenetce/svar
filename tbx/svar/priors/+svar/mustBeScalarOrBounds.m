function mustBeScalarOrBounds(x)
%MUSTBESCALARORBOUNDS Validator for fixed scalars or optimizer bounds.

if isa(x, "hyperprior")
    return
end

if numel(x) ~= 1 && numel(x) ~= 2
    error("mustBeScalarOrBounds:invalidHyperparam", ...
        "Value must be a scalar or a 2-element [lower upper] bound.");
end

if any(x <= 0)
    error("mustBeScalarOrBounds:invalidHyperparam", ...
        "Value must be positive.");
end

if numel(x) == 2
    if any(~isfinite(x))
        error("mustBeScalarOrBounds:invalidHyperparam", ...
            "Bounds [lower upper] must both be finite.");
    end
    if x(1) >= x(2)
        error("mustBeScalarOrBounds:invalidHyperparam", ...
            "Bounds must satisfy lower < upper.");
    end
end
end
