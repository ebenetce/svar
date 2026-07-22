classdef minnesotaPriorClassesTest < matlab.unittest.TestCase
    %minnesotaPriorClassesTest Tests for Minnesota prior class dispatch.

    methods (Test)
        function concreteConstructorsReturnExpectedClasses(testCase)
            psi = [1 4 9];

            inwPrior = svar.minnesotainwbvarm(3, 2, psi);
            normalPrior = svar.minnesotanbvarm(3, 2, psi);

            testCase.verifyClass(inwPrior, "svar.minnesotainwbvarm");
            testCase.verifyClass(normalPrior, "svar.minnesotanbvarm");
        end

        function abstractBaseCannotBeInstantiated(testCase)
            testCase.verifyError(@() svar.minnesotabvarmBase(), "MATLAB:class:abstract");
        end

        function abstractSpecBaseCannotBeInstantiated(testCase)
            testCase.verifyError(@() svar.minnesotaBaseSpec(), "MATLAB:class:abstract");
        end

        function concreteSpecConstructorsAreFactoryOnly(testCase)
            testCase.verifyError(@() svar.minnesotainwSpec(), ...
                "MATLAB:class:MethodRestricted");
            testCase.verifyError(@() svar.minnesotanSpec(), ...
                "MATLAB:class:MethodRestricted");
        end

        function inwCovarianceUsesFullCoefficientLayout(testCase)
            prior = svar.minnesotainwbvarm(3, 2, [1 4 9]);
            numCoefficients = minnesotaPriorClassesTest.numEquationCoefficients(prior);

            testCase.verifySize(prior.Mu, [numCoefficients*prior.NumSeries 1]);
            testCase.verifySize(prior.V, ...
                [numCoefficients*prior.NumSeries numCoefficients*prior.NumSeries]);
            testCase.verifySize(prior.Omega, [prior.NumSeries prior.NumSeries]);
        end

        function normalCovarianceUsesFullCoefficientLayout(testCase)
            prior = svar.minnesotanbvarm(3, 2, [1 4 9]);
            numCoefficients = minnesotaPriorClassesTest.numEquationCoefficients(prior);

            testCase.verifySize(prior.Mu, [numCoefficients*prior.NumSeries 1]);
            testCase.verifySize(prior.V, ...
                [numCoefficients*prior.NumSeries numCoefficients*prior.NumSeries]);
        end

        function normalPriorSetsFixedSigmaFromPpsi(testCase)
            prior = svar.minnesotanbvarm(3, 2, [1 4 9]);

            testCase.verifyEqual(prior.Sigma, diag([1 4 9]), AbsTol=0);
        end

        function normalPriorShrinksOnlyCrossLagCoefficientsWithLambda2(testCase)
            prior = svar.minnesotanbvarm(2, 1, [1 4], ...
                lambda1=0.2, lambda2=0.5, lambda3=1);

            testCase.verifyEqual(diag(prior.V), ...
                [0.04; 0.0025; 1e4; 0.04; 0.04; 1e4], AbsTol=1e-14);
        end

        function normalPriorEstimateReturnsNormalPosterior(testCase)
            prior = svar.minnesotanbvarm(2, 1, [1 4], ...
                lambda1=0.2, lambda2=0.5, lambda3=1);
            Y = minnesotaPriorClassesTest.normalLikelihoodData();

            posterior = estimate(prior, Y, Display="off");

            testCase.verifyClass(posterior, "normalbvarm");
        end

        function normalPriorMarginalLikelihoodMethodsAgree(testCase)
            prior = svar.minnesotanbvarm(2, 1, [1 4], ...
                lambda1=0.2, lambda2=0.5, lambda3=1);
            Y = minnesotaPriorClassesTest.normalLikelihoodData();

            logML = logMarginalLikelihood(prior, Y);
            ml = marginalLikelihood(prior, Y);

            testCase.verifyTrue(isfinite(logML));
            testCase.verifyEqual(ml, exp(logML), AbsTol=1e-12);
        end

        function factoryDispatchesExpectedSpecClasses(testCase)
            inwSpec = minnesotaSpec("inw");
            normalSpec = minnesotaSpec("normal");
            legacySpec = minnesotaSpec("inw", lambda1=[1e-3 5]);

            testCase.verifyClass(inwSpec, "svar.minnesotainwSpec");
            testCase.verifyClass(normalSpec, "svar.minnesotanSpec");
            testCase.verifyEqual(legacySpec.lambda1, [1e-3 5], AbsTol=0);
        end

        function factoryDispatchesAliases(testCase)
            testCase.verifyClass(minnesotaSpec("semiconjugate"), "svar.minnesotainwSpec");
            testCase.verifyClass(minnesotaSpec("fixedsigma"), "svar.minnesotanSpec");
        end

        function retiredConjugateAliasesPointAtTheClass(testCase)
            % The conjugate family left the spec API: minnesotamniwbvarm
            % takes the data directly and glp tunes it without a spec.
            for alias = ["mniw" "conjugate" "matrixnormal"]
                testCase.verifyError(@() minnesotaSpec(alias), ...
                    "minnesotaSpec:retiredMethod");
            end
        end

        function unknownMethodMessageListsAliases(testCase)
            import matlab.unittest.constraints.ContainsSubstring

            message = minnesotaPriorClassesTest.unknownMethodMessage();

            testCase.verifyThat(message, ...
                ContainsSubstring("Available Method aliases:"));
            testCase.verifyThat(message, ...
                ContainsSubstring("inw    (Independent Normal-Wishart):"));
            testCase.verifyThat(message, ...
                ContainsSubstring("normal (fixed-Sigma Normal):"));
            testCase.verifyThat(message, ...
                ContainsSubstring("Hyphens, underscores, and spaces are ignored"));
        end

        function specsInheritSharedBase(testCase)
            testCase.verifyTrue(isa(minnesotaSpec("inw"), "svar.minnesotaBaseSpec"));
            testCase.verifyTrue(isa(minnesotaSpec("normal"), "svar.minnesotaBaseSpec"));
        end

        function specsBuildExpectedModelClasses(testCase)
            psi = [1 4 9];

            inwPrior = minnesotaSpec("inw").build(3, 2, psi);
            normalPrior = minnesotaSpec("normal").build(3, 2, psi);

            testCase.verifyClass(inwPrior, "svar.minnesotainwbvarm");
            testCase.verifyClass(normalPrior, "svar.minnesotanbvarm");
        end

        function inwSpecPacksHyperpriorProperties(testCase)
            spec = minnesotaSpec("inw", ...
                lambda1=hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5]), ...
                lambda2=hyperprior("Gamma", 0.5, 0.25, Bounds=[1e-4 1]), ...
                lambda3=1);

            [x0, lb, ub, names] = spec.pack();

            testCase.verifyEqual(names, ["lambda1" "lambda2"]);
            testCase.verifyEqual(x0, [0.2 0.5], AbsTol=1e-14);
            testCase.verifyEqual(lb, [1e-4 1e-4], AbsTol=0);
            testCase.verifyEqual(ub, [5 1], AbsTol=0);
        end

        function normalSpecPacksHyperpriorProperties(testCase)
            spec = minnesotaSpec("normal", ...
                lambda1=hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5]), ...
                lambda2=hyperprior("Gamma", 0.5, 0.25, Bounds=[1e-4 1]), ...
                lambda3=1);

            [x0, lb, ub, names] = spec.pack();

            testCase.verifyEqual(names, ["lambda1" "lambda2"]);
            testCase.verifyEqual(x0, [0.2 0.5], AbsTol=1e-14);
            testCase.verifyEqual(lb, [1e-4 1e-4], AbsTol=0);
            testCase.verifyEqual(ub, [5 1], AbsTol=0);
        end

        function specLogHyperpriorUsesPackVectorAndNames(testCase)
            spec = minnesotaSpec("inw", ...
                lambda1=hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5]), ...
                lambda2=[1e-4 1], ...
                lambda3=hyperprior("Gamma", 1, 1, Bounds=[1e-4 50]));
            [x0, ~, ~, names] = spec.pack();

            logp = spec.logHyperprior(x0, names);
            expected = spec.lambda1.logpdf(x0(1)) ...
                + spec.lambda3.logpdf(x0(3));

            testCase.verifyEqual(logp, expected, AbsTol=1e-14);
        end
    end

    methods (Static, Access = private)
        function n = numEquationCoefficients(prior)
            n = prior.NumSeries*prior.P ...
                + double(prior.IncludeConstant) ...
                + double(prior.IncludeTrend) ...
                + prior.NumPredictors;
        end

        function message = unknownMethodMessage()
            try
                minnesotaSpec("unknown");
                message = "";
            catch ME
                message = string(ME.message);
            end
        end

        function Y = normalLikelihoodData()
            Y = [ ...
                1.0  2.0
                1.4  2.1
                1.9  2.5
                2.2  2.9
                2.8  3.0
                3.1  3.4
                3.6  3.7
                4.0  4.1];
        end
    end
end
