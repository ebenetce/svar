classdef minnesotaPriorClassesTest < matlab.unittest.TestCase
    %minnesotaPriorClassesTest Tests for Minnesota prior class dispatch.

    methods (TestClassSetup)
        function addProjectToPath(testCase)
            projectFolder = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectFolder, "tbx", "svar")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectFolder, "tbx", "svar", "priors")));
        end
    end

    methods (Test)
        function concreteConstructorsReturnExpectedClasses(testCase)
            psi = [1 4 9];

            mniwPrior = minnesotamniwbvarm(3, 2, ResidualVariances=psi);
            inwPrior = minnesotainwbvarm(3, 2, ResidualVariances=psi);
            normalPrior = minnesotanbvarm(3, 2, ResidualVariances=psi);

            testCase.verifyClass(mniwPrior, "minnesotamniwbvarm");
            testCase.verifyClass(inwPrior, "minnesotainwbvarm");
            testCase.verifyClass(normalPrior, "minnesotanbvarm");
        end

        function abstractBaseCannotBeInstantiated(testCase)
            testCase.verifyError(@() minnesotabvarm(), "MATLAB:class:abstract");
        end

        function inwCovarianceUsesFullCoefficientLayout(testCase)
            prior = minnesotainwbvarm(3, 2, ResidualVariances=[1 4 9]);
            numCoefficients = minnesotaPriorClassesTest.numEquationCoefficients(prior);

            testCase.verifySize(prior.Mu, [numCoefficients*prior.NumSeries 1]);
            testCase.verifySize(prior.V, ...
                [numCoefficients*prior.NumSeries numCoefficients*prior.NumSeries]);
            testCase.verifySize(prior.Omega, [prior.NumSeries prior.NumSeries]);
        end

        function normalCovarianceUsesFullCoefficientLayout(testCase)
            prior = minnesotanbvarm(3, 2, ResidualVariances=[1 4 9]);
            numCoefficients = minnesotaPriorClassesTest.numEquationCoefficients(prior);

            testCase.verifySize(prior.Mu, [numCoefficients*prior.NumSeries 1]);
            testCase.verifySize(prior.V, ...
                [numCoefficients*prior.NumSeries numCoefficients*prior.NumSeries]);
        end

        function normalPriorSetsFixedSigmaFromPpsi(testCase)
            prior = minnesotanbvarm(3, 2, ResidualVariances=[1 4 9]);

            testCase.verifyEqual(prior.Sigma, diag([1 4 9]), AbsTol=0);
        end

        function normalPriorShrinksOnlyCrossLagCoefficientsWithLambda2(testCase)
            prior = minnesotanbvarm(2, 1, ResidualVariances=[1 4], ...
                lambda1=0.2, lambda2=0.5, lambda3=1);

            testCase.verifyEqual(diag(prior.V), ...
                [0.04; 0.0025; 1e4; 0.04; 0.04; 1e4], AbsTol=1e-14);
        end

        function factoryDispatchesExpectedSpecClasses(testCase)
            defaultSpec = minnesotaSpec();
            mniwSpec = minnesotaSpec("mniw");
            inwSpec = minnesotaSpec("inw");
            normalSpec = minnesotaSpec("normal");
            legacySpec = minnesotaSpec(lambda1=[1e-3 5]);

            testCase.verifyClass(defaultSpec, "minnesotamniwSpec");
            testCase.verifyClass(mniwSpec, "minnesotamniwSpec");
            testCase.verifyClass(inwSpec, "minnesotainwSpec");
            testCase.verifyClass(normalSpec, "minnesotanSpec");
            testCase.verifyEqual(legacySpec.lambda1, [1e-3 5], AbsTol=0);
        end

        function specsBuildExpectedModelClasses(testCase)
            psi = [1 4 9];

            mniwPrior = minnesotaSpec().build(3, 2, psi);
            inwPrior = minnesotaSpec("inw").build(3, 2, psi);
            normalPrior = minnesotaSpec("normal").build(3, 2, psi);

            testCase.verifyClass(mniwPrior, "minnesotamniwbvarm");
            testCase.verifyClass(inwPrior, "minnesotainwbvarm");
            testCase.verifyClass(normalPrior, "minnesotanbvarm");
        end
    end

    methods (Static, Access = private)
        function n = numEquationCoefficients(prior)
            n = prior.NumSeries*prior.P ...
                + double(prior.IncludeConstant) ...
                + double(prior.IncludeTrend) ...
                + prior.NumPredictors;
        end
    end
end
