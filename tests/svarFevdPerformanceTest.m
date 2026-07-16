classdef svarFevdPerformanceTest < matlab.perftest.TestCase
    %svarFevdPerformanceTest Benchmarks for FEVD calculations.

    properties (MethodSetupParameter)
        Problem = struct("medium", struct( ...
            "NumSeries", 6, "NumLags", 4, "Horizon", 40))
    end

    properties
        Mdl
        Impact
        Horizon
        ArLagOp
        MaLagOp
    end

    methods (TestMethodSetup)
        function setupProblem(testCase, Problem)
            [testCase.Mdl, testCase.Impact] = ...
                svarFevdPerformanceTest.deterministicVarModel( ...
                Problem.NumSeries, Problem.NumLags);
            testCase.Horizon = Problem.Horizon;
            [testCase.ArLagOp, testCase.MaLagOp] = ...
                svarFevdPerformanceTest.armaLagOperators(testCase.Mdl);
        end
    end

    methods (Test, TestTags = {'Performance'})
        function testSvarFevd(testCase)
            mdl = testCase.Mdl;
            impact = testCase.Impact;
            horizon = testCase.Horizon;
            while testCase.keepMeasuring
                decomposition = svar.fevd(mdl, impact, horizon);
            end
            testCase.verifySize(decomposition, ...
                [horizon + 1 mdl.NumSeries mdl.NumSeries]);
        end

        function testArmaFevdOrthogonalized(testCase)
            arLagOp = testCase.ArLagOp;
            maLagOp = testCase.MaLagOp;
            covariance = testCase.Mdl.Covariance;
            horizon = testCase.Horizon;
            numSeries = testCase.Mdl.NumSeries;
            while testCase.keepMeasuring
                decomposition = armafevd(arLagOp, maLagOp, ...
                    InnovCov=covariance, NumObs=horizon + 1, ...
                    Method="orthogonalized");
            end
            testCase.verifySize(decomposition, ...
                [horizon + 1 numSeries numSeries]);
        end

        function testArmaFevdGeneralized(testCase)
            arLagOp = testCase.ArLagOp;
            maLagOp = testCase.MaLagOp;
            covariance = testCase.Mdl.Covariance;
            horizon = testCase.Horizon;
            numSeries = testCase.Mdl.NumSeries;
            while testCase.keepMeasuring
                decomposition = armafevd(arLagOp, maLagOp, ...
                    InnovCov=covariance, NumObs=horizon + 1, ...
                    Method="generalized");
            end
            testCase.verifySize(decomposition, ...
                [horizon + 1 numSeries numSeries]);
        end

        function testVarmFevdOrthogonalized(testCase)
            mdl = testCase.Mdl;
            horizon = testCase.Horizon;
            while testCase.keepMeasuring
                decomposition = fevd(mdl, NumObs=horizon + 1, ...
                    Method="orthogonalized");
            end
            testCase.verifySize(decomposition, ...
                [horizon + 1 mdl.NumSeries mdl.NumSeries]);
        end

        function testVarmFevdGeneralized(testCase)
            mdl = testCase.Mdl;
            horizon = testCase.Horizon;
            while testCase.keepMeasuring
                decomposition = fevd(mdl, NumObs=horizon + 1, ...
                    Method="generalized");
            end
            testCase.verifySize(decomposition, ...
                [horizon + 1 mdl.NumSeries mdl.NumSeries]);
        end
    end

    methods (Static, Access = private)
        function [mdl, impact] = deterministicVarModel(numSeries, numLags)
            coefficientScale = 0.12/numLags;
            base = reshape(sin(1:numSeries^2), numSeries, numSeries);
            ar = cell(1, numLags);
            for lag = 1:numLags
                ar{lag} = coefficientScale*cos(lag)*base;
                ar{lag} = ar{lag} + (0.18/lag)*eye(numSeries);
            end

            impact = tril(0.15*reshape(cos(1:numSeries^2), ...
                numSeries, numSeries));
            impact = impact + diag(0.8 + 0.1*(1:numSeries));

            mdl = varm(numSeries, numLags);
            mdl.AR = ar;
            mdl.Constant = zeros(numSeries, 1);
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
