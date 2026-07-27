function X = seasonalDummies(numObs, period, nvp)
%seasonalDummies - Create seasonal indicator variables
%   X = seasonalDummies(numObs,PERIOD) creates a seasonal indicator
%   matrix with numObs rows and PERIOD seasonal categories. The first
%   observation is category 1, and categories repeat cyclically.
%
%   X = seasonalDummies(...,DropLast=TF) also specifies whether to omit
%   one category. The default is true, which avoids collinearity with a
%   model intercept.
%
%   X = seasonalDummies(...,Reference=REF) also specifies the category
%   omitted when DropLast is true. REF defaults to PERIOD.
%
%   Example: Monthly dummies with January omitted
%       X = svar.seasonalDummies(36,12);
%
%   See also datetime, month, eye

arguments
    numObs (1,1) {mustBePositive, mustBeInteger}
    period (1,1) {mustBePositive, mustBeInteger}
    nvp.DropLast (1,1) logical = true
    nvp.Reference (1,1) {mustBePositive, mustBeInteger} = period
end

if nvp.Reference > period
    error("seasonalDummies:InvalidReference", ...
        "Reference must be less than or equal to period.");
end

seasons = mod((0:numObs-1)', period) + 1;
indicator = eye(period);
X = indicator(seasons, :);

if nvp.DropLast
    X(:, nvp.Reference) = [];
end
end
