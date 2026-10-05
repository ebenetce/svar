function [impact, structuralMatrix, shockCovariance] = reducedIdentification(covariance)
%REDUCEDIDENTIFICATION Return the unorthogonalized reduced-form innovations.

numSeries = size(covariance, 1);
impact = eye(numSeries);
structuralMatrix = impact;
shockCovariance = covariance;
end
