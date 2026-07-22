classdef glpHyperpriorTest < matlab.unittest.TestCase
    %glpHyperpriorTest Tests for GLP hyperparameter and Psi inputs.

    methods (Test)
        function acceptsFixedNumericPsi(testCase)
            Y = glpHyperpriorTest.sampleData();

            [mdl, info] = glp(2, 1, Y, Psi=[0.5 1.0], ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyClass(mdl, "svar.minnesotamniwbvarm");
            testCase.verifyFalse(info.PsiFree);
            testCase.verifyEqual(info.FinalPsi, [0.5 1.0], AbsTol=0);
        end

        function acceptsNumericPsiBounds(testCase)
            Y = glpHyperpriorTest.sampleData();

            [mdl, info] = glp(2, 1, Y, Psi=[0.1 0.2; 2.0 3.0], ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyClass(mdl, "svar.minnesotamniwbvarm");
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

            testCase.verifyClass(mdl, "svar.minnesotamniwbvarm");
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

        function noChainUnlessNumDrawsRequested(testCase)
            Y = glpHyperpriorTest.persistentData();

            [~, ~, chain] = glp(2, 1, Y, Psi=[0.5 1.0], lambda1=[1e-3 5], ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyEmpty(chain);
        end

        function chainHasOneDrawPerKeptIteration(testCase)
            Y = glpHyperpriorTest.persistentData();
            rng(0, "twister");

            [mdl, ~, chain] = glp(2, 1, Y, Psi=[0.5 1.0], lambda1=[1e-3 5], ...
                lambda4=[1e-3 50], NumDraws=40, BurnIn=10, ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifySize(chain.lambda1, [40 1]);
            testCase.verifySize(chain.Psi, [40 2]);
            testCase.verifySize(chain.Sigma, [2 2 40]);
            testCase.verifySize(chain.Coefficients, ...
                [numel(mdl.Mu)/2, 2, 40]);
            testCase.verifySize(chain.LogPosterior, [40 1]);
            testCase.verifyEqual(chain.BurnIn, 10);
        end

        function chainStaysInsideTheSearchBox(testCase)
            % Candidates outside the box must be rejected, not clipped.
            Y = glpHyperpriorTest.persistentData();
            rng(0, "twister");

            [~, info, chain] = glp(2, 1, Y, Psi=[0.5 1.0], lambda1=[0.05 0.5], ...
                NumDraws=60, BurnIn=10, ProposalScale=3, ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyGreaterThanOrEqual(chain.lambda1, info.LowerBound(1));
            testCase.verifyLessThanOrEqual(chain.lambda1, info.UpperBound(1));
        end

        function fixedHyperparametersAreConstantAlongTheChain(testCase)
            Y = glpHyperpriorTest.persistentData();
            rng(0, "twister");

            [~, ~, chain] = glp(2, 1, Y, Psi=[0.5 1.0], lambda1=[1e-3 5], ...
                lambda3=1, NumDraws=30, BurnIn=5, ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyEqual(chain.lambda3, ones(30, 1), AbsTol=0);
            testCase.verifyEqual(chain.Psi, repmat([0.5 1.0], 30, 1), AbsTol=0);
            testCase.verifyGreaterThan(std(chain.lambda1), 0);
        end

        function samplingWithoutFreeHyperparametersErrors(testCase)
            Y = glpHyperpriorTest.sampleData();

            testCase.verifyError( ...
                @() glp(2, 1, Y, Psi=[0.5 1.0], NumDraws=10), ...
                "glp:nothingToSample");
        end

        function acceptanceRateFallsAsTheProposalWidens(testCase)
            Y = glpHyperpriorTest.persistentData();
            options = glpHyperpriorTest.optimOptions();

            rng(0, "twister");
            [~, ~, tight] = glp(2, 1, Y, Psi=[0.5 1.0], lambda1=[1e-3 5], ...
                NumDraws=200, BurnIn=50, ProposalScale=0.3, OptimOptions=options);
            rng(0, "twister");
            [~, ~, wide] = glp(2, 1, Y, Psi=[0.5 1.0], lambda1=[1e-3 5], ...
                NumDraws=200, BurnIn=50, ProposalScale=4, OptimOptions=options);

            testCase.verifyGreaterThan(tight.AcceptanceRate, wide.AcceptanceRate);
        end

        function chainCentresOnTheMaximiser(testCase)
            % A correct sampler must put the mode inside its own 95% interval.
            Y = glpHyperpriorTest.persistentData();
            rng(0, "twister");

            [mdl, ~, chain] = glp(2, 1, Y, Psi=[0.5 1.0], lambda1=[1e-3 5], ...
                NumDraws=800, BurnIn=400, ProposalScale=0.8, ...
                OptimOptions=glpHyperpriorTest.optimOptions());

            testCase.verifyGreaterThanOrEqual(mdl.lambda1, ...
                quantile(chain.lambda1, 0.025));
            testCase.verifyLessThanOrEqual(mdl.lambda1, ...
                quantile(chain.lambda1, 0.975));
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
