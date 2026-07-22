classdef logMarginalLikelihoodTest < matlab.unittest.TestCase
    %logMarginalLikelihoodTest Tests for analytic BVAR evidence dispatch.

    methods (Test)
        function conjugateModelReturnsFiniteLogEvidence(testCase)
            Y = logMarginalLikelihoodTest.sampleData();
            mdl = conjugatebvarm(2, 1);

            logML = logMarginalLikelihood(mdl, Y);

            testCase.verifyTrue(isscalar(logML));
            testCase.verifyTrue(isfinite(logML));
        end

        function normalModelReturnsFiniteLogEvidence(testCase)
            Y = logMarginalLikelihoodTest.sampleData();
            mdl = normalbvarm(2, 1);

            logML = logMarginalLikelihood(mdl, Y);

            testCase.verifyTrue(isscalar(logML));
            testCase.verifyTrue(isfinite(logML));
        end

        function marginalLikelihoodExponentiatesLogEvidence(testCase)
            Y = logMarginalLikelihoodTest.sampleData();
            mdl = conjugatebvarm(2, 1);

            logML = logMarginalLikelihood(mdl, Y);
            ml = marginalLikelihood(mdl, Y);

            testCase.verifyEqual(ml, exp(logML), RelTol=1e-12);
        end

        function mniwMinnesotaDelegatesToConjugateDispatcher(testCase)
            Y = logMarginalLikelihoodTest.sampleData();
            prior = minnesotamniwbvarm(2, 1, Y, Psi=[1 4]);
            plain = conjugatebvarm(prior.NumSeries, prior.P, ...
                IncludeConstant=prior.IncludeConstant, ...
                IncludeTrend=prior.IncludeTrend, ...
                NumPredictors=prior.NumPredictors, ...
                SeriesNames=prior.SeriesNames, ...
                Mu=prior.Mu, V=prior.V, Omega=prior.Omega, DoF=prior.DoF);

            logML = logMarginalLikelihood(prior, Y);
            expected = logMarginalLikelihood(plain, Y);

            testCase.verifyEqual(logML, expected, AbsTol=1e-12);
        end

        function normalMinnesotaDelegatesToNormalDispatcher(testCase)
            Y = logMarginalLikelihoodTest.sampleData();
            prior = svar.minnesotanbvarm(2, 1, [1 4]);
            plain = normalbvarm(prior.NumSeries, prior.P, ...
                IncludeConstant=prior.IncludeConstant, ...
                IncludeTrend=prior.IncludeTrend, ...
                NumPredictors=prior.NumPredictors, ...
                SeriesNames=prior.SeriesNames, ...
                Mu=prior.Mu, V=prior.V, Sigma=prior.Sigma);

            logML = logMarginalLikelihood(prior, Y);
            expected = logMarginalLikelihood(plain, Y);

            testCase.verifyEqual(logML, expected, AbsTol=1e-12);
        end

        function mniwMinnesotaDummiesReturnFiniteLogEvidence(testCase)
            Y = logMarginalLikelihoodTest.sampleData();
            prior = minnesotamniwbvarm(2, 1, Y, Psi=[1 4], ...
                lambda4=10, lambda5=5);

            logML = logMarginalLikelihood(prior, Y);

            testCase.verifyTrue(isfinite(logML));
        end

        function usesTheSampleStoredOnTheModel(testCase)
            Y = logMarginalLikelihoodTest.sampleData();
            prior = minnesotamniwbvarm(2, 1, Y, Psi=[1 4], lambda4=10);

            testCase.verifyEqual(logMarginalLikelihood(prior), ...
                logMarginalLikelihood(prior, Y), AbsTol=0);
            testCase.verifyEqual(marginalLikelihood(prior), ...
                marginalLikelihood(prior, Y), AbsTol=0);
        end

        function missingSampleErrorsForModelsWithoutOne(testCase)
            testCase.verifyError(@() logMarginalLikelihood(conjugatebvarm(2, 1)), ...
                "logMarginalLikelihood:missingData");
        end

        function semiconjugateModelErrors(testCase)
            Y = logMarginalLikelihoodTest.sampleData();

            testCase.verifyError( ...
                @() logMarginalLikelihood(semiconjugatebvarm(2, 1), Y), ...
                "logMarginalLikelihood:unsupportedModel");
        end

        function minnesotaInwModelErrors(testCase)
            Y = logMarginalLikelihoodTest.sampleData();
            prior = svar.minnesotainwbvarm(2, 1, [1 4]);

            testCase.verifyError(@() logMarginalLikelihood(prior, Y), ...
                "logMarginalLikelihood:unsupportedModel");
        end

        function improperPriorsError(testCase)
            Y = logMarginalLikelihoodTest.sampleData();

            testCase.verifyError(@() logMarginalLikelihood(diffusebvarm(2, 1), Y), ...
                "logMarginalLikelihood:improperPrior");
            testCase.verifyError(@() logMarginalLikelihood(weakbvarm(2, 1), Y), ...
                "logMarginalLikelihood:improperPrior");
            testCase.verifyError(@() logMarginalLikelihood(uniformirbvarm(2, 1), Y), ...
                "logMarginalLikelihood:improperPrior");
        end
    end

    methods (Static, Access = private)
        function Y = sampleData()
            Y = [ ...
                1.0  2.0
                1.4  2.1
                1.9  2.5
                2.2  2.9
                2.8  3.0
                3.1  3.4
                3.6  3.7
                4.0  4.1];
        end
    end
end
