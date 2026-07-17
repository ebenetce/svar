classdef weakbvarmTest < matlab.unittest.TestCase
    %weakbvarmTest Tests for the weak improper reduced-form prior.

    methods (Test)
        function constructorSetsWeakHyperparameters(testCase)
            n = 3;
            p = 2;
            k = n * p + 1;

            prior = weakbvarm(n, p);

            testCase.verifyClass(prior, "weakbvarm");
            testCase.verifyEqual(prior.DoF, 0);
            testCase.verifyEqual(prior.Omega, zeros(n), AbsTol=0);
            testCase.verifyEqual(prior.Mu, zeros(k * n, 1), AbsTol=0);
            testCase.verifyEqual(diag(prior.V), inf(k, 1));
            testCase.verifyEqual(prior.V(~logical(eye(k))), ...
                zeros(k * (k - 1), 1), AbsTol=0);
        end

        function posteriorMatchesWeakNiwFormulas(testCase)
            n = 2;
            p = 2;
            y = weakbvarmTest.deterministicData(60, n);

            prior = weakbvarm(n, p);
            post = estimate(prior, y, Display="off");
            expected = weakbvarmTest.expectedPosterior(y, p);

            testCase.verifyClass(post, "conjugatebvarm");
            testCase.verifyEqual(post.Mu, expected.Mu, AbsTol=1e-10);
            testCase.verifyEqual(post.V, expected.V, AbsTol=1e-10);
            testCase.verifyEqual(post.Omega, expected.Omega, AbsTol=1e-10);
            testCase.verifyEqual(post.DoF, size(y, 1) - p, AbsTol=0);
        end

        function posteriorDiffersFromDiffusePrior(testCase)
            n = 2;
            p = 2;
            rng(7);
            y = randn(80, n);

            weakPost = estimate(weakbvarm(n, p), y, Display="off");
            diffusePost = estimate(diffusebvarm(n, p), y, Display="off");

            testCase.verifyNotEqual(weakPost.DoF, diffusePost.DoF);
        end
    end

    methods (Static, Access = private)
        function y = deterministicData(numObs, numSeries)
            t = (1:numObs)';
            y = [sin(t / 3) + 0.01 * t, cos(t / 5) - 0.02 * t];
            y = y(:, 1:numSeries);
        end

        function expected = expectedPosterior(y, p)
            yResponse = y((p + 1):end, :);
            design = weakbvarmTest.defaultLagDesign(y, p);
            coefficients = design \ yResponse;
            residuals = yResponse - design * coefficients;

            expected.Mu = coefficients(:);
            expected.V = inv(design' * design);
            expected.Omega = residuals' * residuals;
        end

        function design = defaultLagDesign(y, p)
            numRows = size(y, 1) - p;
            numSeries = size(y, 2);
            design = zeros(numRows, numSeries * p + 1);

            for lag = 1:p
                columns = ((lag - 1) * numSeries + 1):(lag * numSeries);
                rows = (p + 1 - lag):(size(y, 1) - lag);
                design(:, columns) = y(rows, :);
            end

            design(:, end) = 1;
        end
    end
end
