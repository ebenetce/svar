function mustBeScalarOrBounds(x)
%MUSTBESCALARORBOUNDS Validator: x is a positive scalar (fixed value, Inf
%   allowed as an "off" marker where the caller's convention uses one - e.g.
%   MINNESOTASPEC.lambda4/lambda5) or a finite 2-element [lower upper] with
%   lower < upper (a free bound for optimisation).
%
%   Shared by MINNESOTASPEC and MINNESOTAINWSPEC so both hyperparameter recipe
%   classes validate the scalar-or-bounds convention identically - one
%   implementation, not a copy per class.
%
%   See also MINNESOTASPEC, MINNESOTAINWSPEC.

if isa(x, 'hyperprior')
    return
end

if numel(x) ~= 1 && numel(x) ~= 2
    error("mustBeScalarOrBounds:invalidHyperparam", ...
        "Value must be a scalar (fixed) or a 2-element [lower upper] bound " + ...
        "(free); got %d elements.", numel(x));
end
if any(x <= 0)
    error("mustBeScalarOrBounds:invalidHyperparam", ...
        "Value must be positive (Inf permitted as a scalar 'off' marker where " + ...
        "the caller's convention defines one).");
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