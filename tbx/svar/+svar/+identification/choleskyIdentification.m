function [impact, structuralMatrix, shockCovariance] = choleskyIdentification(covariance, draw)
%CHOLESKYIDENTIFICATION Identify unit-variance recursive structural shocks.

arguments
    covariance (:,:) double
    draw (1,1) double {mustBeInteger, mustBePositive}
end

[impact, factorStatus] = chol(covariance, "lower");
if factorStatus ~= 0
    error("svar:identify:NonPositiveDefiniteCovariance", ...
        "Covariance draw %d must be positive definite.", draw);
end

structuralMatrix = impact;
shockCovariance = eye(size(covariance));
end
