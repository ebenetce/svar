classdef choleskyIdentificationTest < matlab.unittest.TestCase
    methods (Test)
        function facadeReturnsCholeskyStructuralDraws(testCase)
            post = struct("SeriesNames", ["1" "2"], "P", 1);
            coeff = zeros(4, 3);
            sigma = repmat([1 0.25; 0.25 2], 1, 1, 3);

            [draws, weights, diagnostics] = svar.identify( ...
                post, "cholesky", coeff, sigma, ...
                Horizon=2, FevdHorizon=2, Verbose=false);

            testCase.verifySize(draws.Impact, [2 2 3]);
            testCase.verifySize(draws.IRF, [3 2 2 3]);
            testCase.verifySize(draws.FEVD, [2 2 2 3]);
            testCase.verifyEqual(draws.ShockOrder, ["1" "2"]);
            testCase.verifyEqual(weights, ones(3, 1)/3, AbsTol=1e-12);
            testCase.verifyEqual(diagnostics.Strategy, "cholesky");
            testCase.verifyEqual(diagnostics.NumAcceptedDraws, 3);
        end

        function classSupportsOrderingAndShockSubset(testCase)
            coeff = zeros(4, 2);
            sigma = repmat([1 0.25; 0.25 2], 1, 1, 2);
            id = svar.identification.choleskyIdentification(["1" "2"], ...
                NumLags=1, Ordering=["2" "1"]);
            id.Verbose = false;

            [allDraws, ~, allDiagnostics] = id.identify( ...
                coeff, sigma, Horizon=2, FevdHorizon=2);
            [subsetDraws, ~, subsetDiagnostics] = id.identify( ...
                coeff, sigma, Horizon=2, FevdHorizon=2, Shocks="1");

            testCase.verifyEqual(allDraws.ShockOrder, ["2" "1"]);
            testCase.verifyEqual(allDiagnostics.Ordering, [2 1]);
            testCase.verifyEqual(subsetDraws.ShockOrder, "1");
            testCase.verifyEqual(subsetDiagnostics.SelectedShockIndices, 2);
            testCase.verifyEqual(subsetDraws.Impact(:, 1, :), ...
                allDraws.Impact(:, 2, :), AbsTol=1e-12);
        end

        function facadeRejectsUnsupportedMethod(testCase)
            post = struct("SeriesNames", ["1" "2"], "P", 1);
            coeff = zeros(4, 1);
            sigma = eye(2);

            testCase.verifyError(@() svar.identify( ...
                post, "other", coeff, sigma), ...
                "svar:identify:UnsupportedMethod");
        end
    end
end
