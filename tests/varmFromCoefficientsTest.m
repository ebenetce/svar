classdef varmFromCoefficientsTest < matlab.unittest.TestCase
    %varmFromCoefficientsTest Tests BVAR draw reconstruction as VARM.

    methods (Test)
        function leavesBetaEmptyWithoutPredictors(testCase)
            template = conjugatebvarm(2, 1, IncludeConstant=true, ...
                IncludeTrend=false, SeriesNames=["Y1" "Y2"]);
            coefficients = [ ...
                0.11 0.12
                0.21 0.22
                1.10 1.20];
            covariance = [1.0 0.1; 0.1 0.9];

            mdl = svar.varmFromCoefficients(template, coefficients, covariance);

            testCase.verifyEqual(mdl.AR{1}, coefficients(1:2, :).', AbsTol=0);
            testCase.verifyEqual(mdl.Constant, coefficients(3, :).', AbsTol=0);
            testCase.verifySize(mdl.Beta, [2 0]);
        end

        function mapsPredictorRowsToBeta(testCase)
            template = conjugatebvarm(2, 2, IncludeConstant=true, ...
                IncludeTrend=true, NumPredictors=2, SeriesNames=["Y1" "Y2"]);
            coefficients = [ ...
                0.11 0.12
                0.21 0.22
                0.31 0.32
                0.41 0.42
                1.10 1.20
                2.10 2.20
                3.10 3.20
                4.10 4.20];
            covariance = [1.0 0.2; 0.2 0.8];

            mdl = svar.varmFromCoefficients(template, coefficients, covariance);

            testCase.verifyEqual(mdl.SeriesNames, ["Y1" "Y2"]);
            testCase.verifyEqual(mdl.AR{1}, coefficients(1:2, :).', AbsTol=0);
            testCase.verifyEqual(mdl.AR{2}, coefficients(3:4, :).', AbsTol=0);
            testCase.verifyEqual(mdl.Constant, coefficients(5, :).', AbsTol=0);
            testCase.verifyEqual(mdl.Trend, coefficients(6, :).', AbsTol=0);
            testCase.verifyEqual(mdl.Beta, coefficients(7:8, :).', AbsTol=0);
            testCase.verifyEqual(mdl.Covariance, covariance, AbsTol=0);
        end

        function sizeMismatchIncludesPredictorRows(testCase)
            template = conjugatebvarm(2, 2, NumPredictors=2);
            coefficients = zeros(5, 2);
            covariance = eye(2);

            testCase.verifyError( ...
                @() svar.varmFromCoefficients(template, coefficients, covariance), ...
                "varmFromCoefficients:sizeMismatch");
        end
    end
end
