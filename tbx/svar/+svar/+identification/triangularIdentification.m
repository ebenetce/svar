function [impact, structuralMatrix, shockCovariance] = triangularIdentification(covariance, draw)
%TRIANGULARIDENTIFICATION Identify BEAR-style triangular structural shocks.

[impact, ~, ~] = svar.identification.choleskyIdentification(covariance, draw);
diagonal = diag(impact);
structuralMatrix = impact./diagonal.';
shockCovariance = diag(diagonal.^2);
end
