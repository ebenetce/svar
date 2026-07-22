classdef minnesotamniwbvarmTest < matlab.unittest.TestCase
    %minnesotamniwbvarmTest Tests for the conjugate Minnesota prior class.

    methods (Test)

        % ---- construction ------------------------------------------------

        function estimatesPsiFromDataByDefault(testCase)
            Y = minnesotamniwbvarmTest.sampleData();

            prior = svar.minnesotamniwbvarm(2, 1, Y);

            testCase.verifyEqual(prior.ResidualVariances, ...
                estimateResidualVariances(Y, 1, Method="exact"), AbsTol=1e-12);
        end

        function acceptsNamedPsiEstimator(testCase)
            Y = minnesotamniwbvarmTest.sampleData();

            prior = svar.minnesotamniwbvarm(2, 1, Y, Psi="conditional");

            testCase.verifyEqual(prior.ResidualVariances, ...
                estimateResidualVariances(Y, 1, Method="conditional"), AbsTol=1e-12);
        end

        function acceptsNumericPsi(testCase)
            Y = minnesotamniwbvarmTest.sampleData();

            prior = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4]);

            testCase.verifyEqual(prior.ResidualVariances, [1 4], AbsTol=0);
            testCase.verifyEqual(prior.Omega, diag([1 4]), AbsTol=0);
        end

        function storesSampleAndAdoptsTableNames(testCase)
            Y = minnesotamniwbvarmTest.sampleData();
            T = array2table(Y, VariableNames=["GDP" "CPI"]);

            prior = svar.minnesotamniwbvarm(2, 1, T, Psi=[1 4]);

            testCase.verifyEqual(string(prior.SeriesNames), ["GDP" "CPI"]);
            testCase.verifyEqual(prior.Y, Y, AbsTol=0);
        end

        function rejectsInvalidSample(testCase)
            Y = minnesotamniwbvarmTest.sampleData();

            testCase.verifyError(@() svar.minnesotamniwbvarm(3, 1, Y), ...
                "minnesotamniwbvarm:invalidData");
            testCase.verifyError(@() svar.minnesotamniwbvarm(2, 1, Y(1,:)), ...
                "minnesotamniwbvarm:invalidData");
        end

        function rejectsInvalidPsi(testCase)
            Y = minnesotamniwbvarmTest.sampleData();

            testCase.verifyError(@() svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 2 3]), ...
                "minnesotamniwbvarm:invalidPsi");
            testCase.verifyError(@() svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 -1]), ...
                "minnesotamniwbvarm:invalidPsi");
        end

        % ---- the prior is fully materialised at construction --------------

        function dummiesAreFoldedInByTheConstructor(testCase)
            Y = minnesotamniwbvarmTest.sampleData();

            off = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4]);
            on  = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4], lambda4=1, lambda5=2);

            testCase.verifyEqual(off.NumDummyObservations, 0);
            testCase.verifyEqual(on.NumDummyObservations, 3);   % n rows + 1 row
            testCase.verifyEqual(on.DoF, off.DoF + 3, AbsTol=0);
            testCase.verifyNotEqual(on.V, off.V);

            % The prior MEAN is untouched at PriorMean = 1: both dummies
            % encode the same unit-root / co-persistence belief the mean
            % already states, so they are satisfied exactly at B0 and only
            % sharpen V and DoF. Move the mean off 1 and they do shift it.
            testCase.verifyEqual(on.Mu, off.Mu, AbsTol=0);
        end

        function dummiesShiftTheMeanAwayFromAUnitRootPriorMean(testCase)
            Y = minnesotamniwbvarmTest.sampleData();

            off = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4], PriorMean=0.5);
            on  = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4], PriorMean=0.5, ...
                lambda4=1, lambda5=2);

            testCase.verifyNotEqual(on.Mu, off.Mu);
        end

        function marginalLikelihoodSeesTheDummies(testCase)
            % Regression guard. The dummy observations previously reached the
            % marginal likelihood through a separate late-applied code path,
            % so a consumer that missed it silently scored lambda4/lambda5 as
            % if they were off.
            Y = minnesotamniwbvarmTest.sampleData();

            off = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4]);
            on  = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4], lambda4=1);

            testCase.verifyNotEqual( ...
                logMarginalLikelihood(on), logMarginalLikelihood(off));
        end

        function lambdasOffByDefault(testCase)
            Y = minnesotamniwbvarmTest.sampleData();

            prior = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4]);

            testCase.verifyEqual(prior.lambda4, Inf);
            testCase.verifyEqual(prior.lambda5, Inf);
            testCase.verifyEqual(prior.lambda2, 1);   % pinned under conjugacy
        end

        function priorMomentsMatchMinnesotaTarget(testCase)
            Y = minnesotamniwbvarmTest.sampleData();
            psi = [1 4];

            prior = svar.minnesotamniwbvarm(2, 2, Y, Psi=psi, ...
                lambda1=0.2, lambda3=1, Vc=1e4);

            % V(lag l, regressor j) = lambda1^2 / (l^(2*lambda3) * psi_j),
            % with the deterministic block at Vc.
            expected = [0.04./psi, 0.04/4./psi, 1e4]';
            testCase.verifyEqual(diag(prior.V), expected, AbsTol=1e-14);
            testCase.verifyEqual(prior.DoF, 4, AbsTol=0);   % n + 2
        end

        % ---- the stored sample drives the inherited surface ---------------

        function estimateWorksWithoutASampleArgument(testCase)
            Y = minnesotamniwbvarmTest.sampleData();
            prior = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4], lambda4=1);

            posterior = estimate(prior);

            testCase.verifyClass(posterior, "conjugatebvarm");
            testCase.verifyEqual(posterior.Mu, ...
                estimate(prior, Y, Display="off").Mu, AbsTol=0);
        end

        function simulateAndForecastWorkWithoutASampleArgument(testCase)
            rng(0);
            Y = minnesotamniwbvarmTest.sampleData();
            prior = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4], lambda4=1);

            [coefficients, covariances] = simulate(prior, NumDraws=5);
            forecasts = forecast(prior, 4);

            testCase.verifySize(coefficients, [6 5]);
            testCase.verifySize(covariances, [2 2 5]);
            testCase.verifySize(forecasts, [4 2]);
        end

        function rejectsASampleOtherThanTheStoredOne(testCase)
            Y = minnesotamniwbvarmTest.sampleData();
            prior = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4]);

            testCase.verifyError(@() estimate(prior, Y(1:10,:)), ...
                "minnesotamniwbvarm:sampleMismatch");
            testCase.verifyError(@() simulate(prior, Y(1:10,:), NumDraws=2), ...
                "minnesotamniwbvarm:sampleMismatch");
        end

        function acceptsTheStoredSampleExplicitly(testCase)
            Y = minnesotamniwbvarmTest.sampleData();
            prior = svar.minnesotamniwbvarm(2, 1, Y, Psi=[1 4]);

            testCase.verifyClass(estimate(prior, Y, Display="off"), ...
                "conjugatebvarm");
        end

    end

    methods (Static, Access = private)
        function Y = sampleData()
            rng(0);
            Y = cumsum(randn(60, 2))/10 + 1;
        end
    end
end
