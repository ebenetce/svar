classdef choleskyIdentificationTest < matlab.unittest.TestCase
    methods (Test)
        function identifiesModelCovariance(testCase)
            mdl = choleskyIdentificationTest.model();

            [impact, info] = svar.identify(mdl, "cholesky");

            testCase.verifyEqual(impact, chol(mdl.Covariance, "lower"), ...
                AbsTol=1e-12);
            testCase.verifyEqual(info.Method, "cholesky");
            testCase.verifyEqual(info.Ordering, [1 2]);
            testCase.verifyEqual(info.ShockOrder, ["Output" "Prices"]);
            testCase.verifyEqual(info.NumDraws, 1);
        end

        function identifiesCovarianceDraws(testCase)
            mdl = choleskyIdentificationTest.model();
            sigma = cat(3, mdl.Covariance, 1.5*mdl.Covariance);

            [impact, info] = svar.identify(mdl, "cholesky", Sigma=sigma);

            testCase.verifySize(impact, [2 2 2]);
            testCase.verifyEqual(info.NumDraws, 2);
            for draw = 1:2
                testCase.verifyEqual(impact(:, :, draw)*impact(:, :, draw).', ...
                    sigma(:, :, draw), AbsTol=1e-12);
            end
        end

        function supportsNamedAndNumericOrderings(testCase)
            mdl = choleskyIdentificationTest.model();
            ordering = [2 1];

            [fromIndices, indexInfo] = svar.identify(mdl, "recursive", ...
                Ordering=ordering);
            [fromNames, nameInfo] = svar.identify(mdl, "triangular", ...
                Ordering=["Prices" "Output"]);
            expected = zeros(2);
            expected(ordering, :) = chol(mdl.Covariance(ordering, ordering), ...
                "lower");

            testCase.verifyEqual(fromIndices, expected, AbsTol=1e-12);
            testCase.verifyEqual(fromNames, expected, AbsTol=1e-12);
            testCase.verifyEqual(indexInfo.Method, "cholesky");
            testCase.verifyEqual(nameInfo.Method, "triangular");
            testCase.verifyEqual(nameInfo.ShockOrder, ["Prices" "Output"]);
            testCase.verifyEqual(fromNames*fromNames.', mdl.Covariance, ...
                AbsTol=1e-12);
        end

        function requiresCovarianceDrawForBayesianVar(testCase)
            Y = randn(40, 2);
            prior = minnesotabvarm(2, 1, Y, SeriesNames=["Output" "Prices"]);
            post = estimate(prior, Display="off");

            testCase.verifyError(@() svar.identify(post, "cholesky"), ...
                "svar:identify:MissingCovarianceDraw");
        end

        function identifiesBayesianCovarianceDraws(testCase)
            Y = randn(40, 2);
            prior = minnesotabvarm(2, 1, Y, SeriesNames=["Output" "Prices"]);
            post = estimate(prior, Display="off");
            [~, sigma] = simulate(post, NumDraws=2);

            [impact, info] = svar.identify(post, "cholesky", Sigma=sigma);

            testCase.verifySize(impact, [2 2 2]);
            testCase.verifyEqual(info.NumDraws, 2);
            for draw = 1:2
                testCase.verifyEqual(impact(:, :, draw)*impact(:, :, draw).', ...
                    sigma(:, :, draw), AbsTol=1e-12);
            end
        end

        function rejectsUnsupportedMethod(testCase)
            mdl = choleskyIdentificationTest.model();

            testCase.verifyError(@() svar.identify(mdl, "sign"), ...
                "svar:identify:UnsupportedMethod");
        end

        function rejectsNonModelInput(testCase)
            testCase.verifyError(@() svar.identify(struct(), "cholesky"), ...
                "MATLAB:validators:mustBeA");
        end

        function rejectsInvalidOrdering(testCase)
            mdl = choleskyIdentificationTest.model();

            testCase.verifyError(@() svar.identify(mdl, "cholesky", ...
                Ordering=[1 1]), "svar:identify:InvalidOrdering");
        end

        function rejectsInvalidCovariance(testCase)
            mdl = choleskyIdentificationTest.model();

            testCase.verifyError(@() svar.identify(mdl, "cholesky", ...
                Sigma=[1 0; 0 0]), ...
                "svar:identify:NonPositiveDefiniteCovariance");
            testCase.verifyError(@() svar.identify(mdl, "cholesky", ...
                Sigma=eye(3)), "svar:identify:InvalidCovariance");
        end
    end

    methods (Static, Access = private)
        function mdl = model()
            mdl = varm(2, 1);
            mdl.SeriesNames = ["Output" "Prices"];
            mdl.Covariance = [1 0.25; 0.25 2];
        end
    end
end
