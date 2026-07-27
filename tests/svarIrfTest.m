classdef svarIrfTest < matlab.unittest.TestCase
    %svarIrfTest Tests for structural VAR impulse responses.

    methods (Test)
        function returnsExpectedShapeAndImpactResponse(testCase)
            [mdl, impact] = svarIrfTest.stableVarModel();
            horizon = 6;

            actual = svar.irf(mdl, impact, horizon);

            testCase.verifySize(actual, [horizon + 1 mdl.NumSeries size(impact, 2)]);
            testCase.verifyEqual(squeeze(actual(1, :, :)), impact, AbsTol=0);
        end

        function reusesPrecomputedCompanionPowers(testCase)
            [mdl, impact] = svarIrfTest.stableVarModel();
            horizon = 6;
            Phi = svar.companionPower(mdl, horizon);
            expected = svar.irf(mdl, impact, horizon);

            actual = svar.irf(mdl, impact, horizon, Phi=Phi);

            testCase.verifyEqual(actual, expected, AbsTol=0);
        end

        function returnsCompanionPowersForReuse(testCase)
            [mdl, impact] = svarIrfTest.stableVarModel();
            horizon = 6;

            [responses, Phi] = svar.irf(mdl, impact, horizon);
            reused = svar.irf(mdl, impact(:, 1), horizon, Phi=Phi);

            testCase.verifyEqual(Phi, svar.companionPower(mdl, horizon), ...
                AbsTol=0);
            testCase.verifyEqual(reused, responses(:, :, 1), AbsTol=0);
        end

        function matchesVarmOrthogonalizedIrf(testCase)
            [mdl, impact] = svarIrfTest.stableVarModel();
            horizon = 6;
            expected = permute(irf(mdl, NumObs=horizon + 1, ...
                Method="orthogonalized"), [1 3 2]);

            actual = svar.irf(mdl, impact, horizon);

            testCase.verifyEqual(actual, expected, AbsTol=1e-12);
        end

        function matchesArmaIrf(testCase)
            [mdl, impact] = svarIrfTest.stableVarModel();
            horizon = 6;
            [arLagOp, maLagOp] = svarIrfTest.armaLagOperators(mdl);
            expected = permute(armairf(arLagOp, maLagOp, ...
                InnovCov=mdl.Covariance, NumObs=horizon + 1, ...
                Method="orthogonalized"), [1 3 2]);

            actual = svar.irf(mdl, impact, horizon);

            testCase.verifyEqual(actual, expected, AbsTol=1e-12);
        end
    end

    methods (Static, Access = private)
        function [mdl, impact] = stableVarModel()
            ar1 = [0.35 0.05 -0.02; 0.01 0.30 0.04; -0.03 0.02 0.25];
            ar2 = [0.08 0.01 0.00; 0.02 0.05 -0.01; 0.00 0.01 0.04];
            impact = [1.0 0.0 0.0; 0.2 0.8 0.0; -0.1 0.1 0.7];

            mdl = varm(3, 2);
            mdl.AR = {ar1, ar2};
            mdl.Constant = zeros(3, 1);
            mdl.Covariance = impact*impact.';
        end

        function [arLagOp, maLagOp] = armaLagOperators(mdl)
            arCoefficients = cell(1, mdl.P + 1);
            arCoefficients{1} = eye(mdl.NumSeries);
            for lag = 1:mdl.P
                arCoefficients{lag + 1} = -mdl.AR{lag};
            end
            arLagOp = LagOp(arCoefficients, Lags=0:mdl.P);
            maLagOp = LagOp(eye(mdl.NumSeries));
        end
    end
end
