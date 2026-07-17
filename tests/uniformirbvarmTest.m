classdef uniformirbvarmTest < matlab.unittest.TestCase
    %uniformirbvarmTest Tests for the uniform IR reduced-form prior.

    methods (Test)
        function constructorSetsImproperPriorHyperparameters(testCase)
            n = 3;
            p = 2;
            k = n*p + 1;

            prior = uniformirbvarm(n, p);

            testCase.verifyClass(prior, "uniformirbvarm");
            testCase.verifyEqual(prior.DeterminantShift, -3);
            testCase.verifyEqual(prior.NumEquationCoefficients, k);
            testCase.verifyEqual(prior.LogDetExponent, 2);
            testCase.verifyEqual(prior.DoF, 2 - 2*k - n);
            testCase.verifyEqual(prior.Mu, zeros(n*k, 1), AbsTol=0);
            testCase.verifyEqual(prior.Omega, zeros(n), AbsTol=0);
            testCase.verifyEqual(diag(prior.V), inf(k, 1));
            testCase.verifyEqual(prior.V(~logical(eye(k))), ...
                zeros(k*(k - 1), 1), AbsTol=0);
        end

        function constructorHonorsModelShapeOptions(testCase)
            n = 2;
            p = 3;
            numPredictors = 2;
            k = n*p + 1 + numPredictors;

            prior = uniformirbvarm(n, p, IncludeConstant=false, ...
                IncludeTrend=true, NumPredictors=numPredictors);

            testCase.verifyEqual(prior.NumEquationCoefficients, k);
            testCase.verifyEqual(prior.LogDetExponent, (k - 3)/2);
            testCase.verifyEqual(prior.DoF, 2 - 2*k - n);
            testCase.verifySize(prior.Mu, [n*k 1]);
            testCase.verifySize(prior.V, [k k]);
        end

        function constructorHonorsDeterminantShift(testCase)
            n = 3;
            p = 1;
            determinantShift = -1;
            k = n*p + 1;

            prior = uniformirbvarm(n, p, ...
                DeterminantShift=determinantShift);

            testCase.verifyEqual(prior.DeterminantShift, determinantShift);
            testCase.verifyEqual(prior.LogDetExponent, ...
                (k + determinantShift)/2);
            testCase.verifyEqual(prior.DoF, ...
                -2*k - n - 1 - determinantShift);
        end

        function posteriorDoFUsesEffectiveSampleSize(testCase)
            n = 2;
            p = 1;
            numObservations = 40;
            y = uniformirbvarmTest.deterministicData(numObservations, n);
            prior = uniformirbvarm(n, p);

            posterior = estimate(prior, y, Display="off");

            expectedDoF = numObservations - p + prior.DoF;
            testCase.verifyTrue(isa(posterior, "conjugatebvarm"));
            testCase.verifyEqual(posterior.DoF, expectedDoF);
        end

        function posteriorMatchesOlsNiwMoments(testCase)
            n = 2;
            p = 2;
            numObservations = 48;
            y = uniformirbvarmTest.deterministicData(numObservations, n);
            prior = uniformirbvarm(n, p);

            posterior = estimate(prior, y, Display="off");
            expected = uniformirbvarmTest.expectedDefaultPosterior(y, p);

            testCase.verifyEqual(posterior.Mu, expected.Mu, AbsTol=1e-10);
            testCase.verifyEqual(posterior.V, expected.V, AbsTol=1e-10);
            testCase.verifyEqual(posterior.Omega, expected.Omega, ...
                AbsTol=1e-10);
        end

        function posteriorDoFWithPresampleUsesYRows(testCase)
            n = 2;
            p = 2;
            numObservations = 35;
            yAll = uniformirbvarmTest.deterministicData(numObservations + p, n);
            y0 = yAll(1:p, :);
            y = yAll((p + 1):end, :);
            prior = uniformirbvarm(n, p);

            posterior = estimate(prior, y, Y0=y0, Display="off");

            expectedDoF = numObservations + prior.DoF;
            testCase.verifyEqual(posterior.DoF, expectedDoF);
        end

        function defaultPaperQuantitiesMatchCorollary(testCase)
            n = 4;
            k = 17;
            p = 4;

            prior = uniformirbvarm(n, p);

            testCase.verifyEqual(prior.NumEquationCoefficients, k);
            testCase.verifyEqual(prior.DoF, -36);
            testCase.verifyEqual(prior.LogDetExponent, 7);
        end

        function customDeterminantShiftChangesPaperQuantities(testCase)
            n = 4;
            k = 17;
            p = 4;
            determinantShift = -1;

            prior = uniformirbvarm(n, p, ...
                DeterminantShift=determinantShift);

            testCase.verifyEqual(prior.NumEquationCoefficients, k);
            testCase.verifyEqual(prior.DoF, -38);
            testCase.verifyEqual(prior.LogDetExponent, 8);
        end
    end

    methods (Static, Access = private)
        function y = deterministicData(numObservations, numSeries)
            t = (1:numObservations)';
            y = [sin(t/3) + 0.01*t, cos(t/5) - 0.02*t];
            y = y(:, 1:numSeries);
        end

        function expected = expectedDefaultPosterior(y, p)
            yResponse = y((p + 1):end, :);
            design = uniformirbvarmTest.defaultLagDesign(y, p);
            coefficients = design\yResponse;
            residuals = yResponse - design*coefficients;

            expected.Mu = coefficients(:);
            expected.V = inv(design'*design);
            expected.Omega = residuals'*residuals;
        end

        function design = defaultLagDesign(y, p)
            numEffectiveObservations = size(y, 1) - p;
            numSeries = size(y, 2);
            design = zeros(numEffectiveObservations, numSeries*p + 1);

            for lag = 1:p
                columns = ((lag - 1)*numSeries + 1):(lag*numSeries);
                rows = (p + 1 - lag):(size(y, 1) - lag);
                design(:, columns) = y(rows, :);
            end

            design(:, end) = 1;
        end
    end
end
