classdef glpHyperpriorTest < matlab.unittest.TestCase
    %glpHyperpriorTest Tests for GLP hyperparameter and Psi inputs.

    methods (Test)
        function acceptsFixedNumericPsi(testCase)
            Y = glpHyperpriorTest.sampleData();

            [mdl, info] = glp(2, 1, Y, Psi=[0.5 1.0], ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyClass(mdl, "minnesotamniwbvarm");
            testCase.verifyFalse(info.PsiFree);
            testCase.verifyEqual(info.FinalPsi, [0.5 1.0], AbsTol=0);
        end

        function acceptsNumericPsiBounds(testCase)
            Y = glpHyperpriorTest.sampleData();

            [mdl, info] = glp(2, 1, Y, Psi=[0.1 0.2; 2.0 3.0], ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyClass(mdl, "minnesotamniwbvarm");
            testCase.verifyTrue(info.PsiFree);
            testCase.verifyEqual(info.PsiNames, ["Psi1" "Psi2"]);
        end

        function acceptsPsiHyperpriorArray(testCase)
            Y = glpHyperpriorTest.sampleData();
            Psi = arrayfun(@(x) hyperprior("InverseGamma", x, 0.5*x, ...
                X0=x, Bounds=[0.1 10]*x), [0.5 1.0]);

            [mdl, info] = glp(2, 1, Y, Psi=Psi, ...
                lambda1=hyperprior("Gamma", 0.2, 0.4, Bounds=[0.05 0.5]), ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyClass(mdl, "minnesotamniwbvarm");
            testCase.verifyTrue(info.PsiFree);
            testCase.verifyTrue(info.UsedHyperprior);
            testCase.verifyEqual(numel(info.X0), 3);
        end

        function resolvesPsiEstimatorNameOnce(testCase)
            Y = glpHyperpriorTest.sampleData();

            [mdl, info] = glp(2, 1, Y, Psi="conditional", ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            % A named estimator is resolved before the search, so Psi is
            % fixed and contributes nothing to the optimiser's box.
            testCase.verifyFalse(info.PsiFree);
            testCase.verifyEqual(info.FinalPsi, mdl.ResidualVariances, AbsTol=0);
            testCase.verifyEqual(info.FinalPsi, ...
                estimateResidualVariances(Y, 1, Method="conditional"), AbsTol=1e-12);
        end

        function scalarHyperparametersAreFixed(testCase)
            Y = glpHyperpriorTest.sampleData();

            [mdl, info] = glp(2, 1, Y, Psi=[0.5 1.0], lambda1=0.33, ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyEmpty(info.FreeLambdas);
            testCase.verifyEmpty(info.X0);
            testCase.verifyEqual(mdl.lambda1, 0.33, AbsTol=0);
        end

        function boundsMakeHyperparameterFree(testCase)
            Y = glpHyperpriorTest.sampleData();

            [~, info] = glp(2, 1, Y, Psi=[0.5 1.0], ...
                lambda1=[0.05 0.5], lambda4=[1e-4 50], ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyEqual(info.FreeLambdas, ["lambda1" "lambda4"]);
            testCase.verifyEqual(info.LowerBound, [0.05 1e-4], AbsTol=0);
            testCase.verifyEqual(info.UpperBound, [0.5 50], AbsTol=0);
        end

        function tunesLambda4ThroughMarginalLikelihood(testCase)
            % Regression guard: lambda4 enters the objective only if the
            % dummy observations reach the marginal likelihood. When they do
            % not, the objective is flat in lambda4 and the optimiser returns
            % the starting value untouched.
            Y = glpHyperpriorTest.persistentData();

            [mdl, info] = glp(2, 1, Y, Psi=[0.5 1.0], lambda4=[1e-3 50], ...
                OptimOptions=optimoptions("fmincon", Display="off"));

            testCase.verifyEqual(info.FreeLambdas, "lambda4");
            testCase.verifyNotEqual(mdl.lambda4, info.X0(1));
            testCase.verifyEqual(mdl.NumDummyObservations, 2);
        end

        function rejectsMalformedHyperparameter(testCase)
            Y = glpHyperpriorTest.sampleData();

            testCase.verifyError(@() glp(2, 1, Y, lambda1=[5 1]), ...
                "glp:invalidHyperparameter");
            testCase.verifyError(@() glp(2, 1, Y, lambda1=[1 2 3]), ...
                "glp:invalidHyperparameter");
            testCase.verifyError(@() glp(2, 1, Y, lambda1=-1), ...
                "glp:invalidHyperparameter");
        end
    end

    methods (Static, Access = private)
        function Y = sampleData()
            Y = [ ...
                1.00  2.00
                1.10  1.95
                1.21  1.90
                1.30  1.86
                1.42  1.80
                1.55  1.75];
        end

        function Y = persistentData()
            rng(0);
            Y = cumsum(randn(80, 2))/10 + 1;
        end

        function opts = optimOptions()
            opts = optimoptions("fmincon", ...
                Display="off", ...
                MaxIterations=1, ...
                MaxFunctionEvaluations=30);
        end
    end
end
