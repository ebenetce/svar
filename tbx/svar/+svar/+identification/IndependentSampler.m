classdef IndependentSampler
    %INDEPENDENTSAMPLER - Draw independently restricted structural models

    methods (Static)
        function [impact, info] = identify(Mdl, restrictions, nvp)
            arguments
                Mdl (1,1) {mustBeA(Mdl, ["varm", "bvar"])}
                restrictions (1,1) svar.identification.RestrictionSet
                nvp.NumDraws (1,1) double {mustBeInteger, mustBePositive}
                nvp.BatchSize (1,1) double {mustBeInteger, mustBePositive}
                nvp.MaxCandidates (1,1) double {mustBeInteger, mustBePositive}
            end

            numSeries = Mdl.NumSeries;
            horizon = restrictions.requiredHorizon();
            numDraws = nvp.NumDraws;
            impact = zeros(numSeries, numSeries, numDraws);
            rotation = zeros(numSeries, numSeries, numDraws);
            sigma = zeros(numSeries, numSeries, numDraws);

            if isa(Mdl, "bvar")
                coefficients = zeros(numel(Mdl.Mu), numDraws);
            else
                coefficients = [];
            end

            numAccepted = 0;
            numCandidates = 0;
            while numAccepted < numDraws && numCandidates < nvp.MaxCandidates
                batchSize = min(nvp.BatchSize, nvp.MaxCandidates - numCandidates);
                [coefficientBatch, sigmaBatch] = localDrawReducedForm(Mdl, batchSize);

                for candidate = 1:batchSize
                    numCandidates = numCandidates + 1;
                    covariance = sigmaBatch(:, :, candidate);
                    [factor, status] = chol(covariance, "lower");
                    if status ~= 0
                        error("svar:identify:NonPositiveDefiniteCovariance", ...
                            "Candidate covariance %d must be positive definite.", ...
                            numCandidates);
                    end

                    [varMdl, coefficient] = localModelForCandidate( ...
                        Mdl, coefficientBatch, covariance, candidate);
                    candidateRotation = svar.identification.restrictedRotation( ...
                        factor, restrictions);
                    candidateImpact = factor*candidateRotation;
                    state = localState(varMdl, candidateImpact, candidateRotation, ...
                        horizon, restrictions);
                    if ~restrictions.isSatisfied(state)
                        continue
                    end

                    numAccepted = numAccepted + 1;
                    impact(:, :, numAccepted) = candidateImpact;
                    rotation(:, :, numAccepted) = candidateRotation;
                    sigma(:, :, numAccepted) = covariance;
                    if isa(Mdl, "bvar")
                        coefficients(:, numAccepted) = coefficient;
                    end

                    if numAccepted == numDraws
                        break
                    end
                end
            end

            if numAccepted < numDraws
                error("svar:identify:MaxCandidatesExceeded", ...
                    "Only %d of %d requested draws were accepted from %d candidates.", ...
                    numAccepted, numDraws, numCandidates);
            end

            info = struct( ...
                "Method", "independent", ...
                "ShockNames", restrictions.ShockOrder, ...
                "NumDraws", numDraws, ...
                "NumCandidates", numCandidates, ...
                "AcceptanceRate", numDraws/numCandidates, ...
                "Restrictions", restrictions, ...
                "Coefficients", coefficients, ...
                "Sigma", sigma, ...
                "Rotation", rotation, ...
                "Impact", impact);
        end
    end
end

function [coefficients, sigma] = localDrawReducedForm(Mdl, batchSize)
if isa(Mdl, "bvar")
    [coefficients, sigma] = simulate(Mdl, NumDraws=batchSize);
else
    coefficients = [];
    sigma = repmat(Mdl.Covariance, 1, 1, batchSize);
end
end

function [varMdl, coefficient] = localModelForCandidate( ...
        Mdl, coefficientBatch, covariance, candidate)
if isa(Mdl, "bvar")
    coefficient = coefficientBatch(:, candidate);
    varMdl = svar.varmFromCoefficients(Mdl, ...
        reshape(coefficient, [], Mdl.NumSeries), covariance);
else
    coefficient = [];
    varMdl = Mdl;
end
end

function state = localState(varMdl, impact, rotation, horizon, restrictions)
state = struct( ...
    "VarMdl", varMdl, ...
    "Impact", impact, ...
    "Rotation", rotation, ...
    "IRF", svar.irf(varMdl, impact, horizon), ...
    "SeriesNames", restrictions.SeriesNames, ...
    "ShockNames", restrictions.ShockOrder);
end
