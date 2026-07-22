classdef minnesotaPriorClassesTest < matlab.unittest.TestCase
    %minnesotaPriorClassesTest Tests for Minnesota prior dispatch and shared base.

    methods (Test)

        % ---- the minnesotabvarm front door -------------------------------

        function defaultTypeIsConjugate(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            prior = minnesotabvarm(2, 1, Y, Psi=[1 4]);

            testCase.verifyClass(prior, "svar.minnesotamniwbvarm");
        end

        function dispatchesEachType(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            testCase.verifyClass(minnesotabvarm(2, 1, Y, Type="mniw", Psi=[1 4]), ...
                "svar.minnesotamniwbvarm");
            testCase.verifyClass(minnesotabvarm(2, 1, Y, Type="inw", Psi=[1 4]), ...
                "svar.minnesotainwbvarm");
            testCase.verifyClass(minnesotabvarm(2, 1, Y, Type="normal", Psi=[1 4]), ...
                "svar.minnesotanbvarm");
        end

        function dispatchesAliasesIgnoringPunctuation(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            testCase.verifyClass(minnesotabvarm(2, 1, Y, Type="conjugate", Psi=[1 4]), ...
                "svar.minnesotamniwbvarm");
            testCase.verifyClass(minnesotabvarm(2, 1, Y, Type="semiconjugate", Psi=[1 4]), ...
                "svar.minnesotainwbvarm");
            testCase.verifyClass(minnesotabvarm(2, 1, Y, Type="fixed-sigma", Psi=[1 4]), ...
                "svar.minnesotanbvarm");
            testCase.verifyClass(minnesotabvarm(2, 1, Y, Type="Kadiyala_Karlsson", Psi=[1 4]), ...
                "svar.minnesotainwbvarm");
        end

        function forwardsFamilySpecificOptions(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            conjugate = minnesotabvarm(2, 1, Y, Psi=[1 4], lambda4=1, lambda5=2);
            independent = minnesotabvarm(2, 1, Y, Type="inw", Psi=[1 4], lambda2=0.25);

            testCase.verifyEqual(conjugate.lambda4, 1);
            testCase.verifyEqual(conjugate.lambda5, 2);
            testCase.verifyEqual(independent.lambda2, 0.25);
        end

        function rejectsOptionsForeignToTheChosenType(testCase)
            % lambda2 is pinned to 1 by conjugacy, so passing it to "mniw"
            % must fail rather than be silently dropped.
            Y = minnesotaPriorClassesTest.sampleData();

            testCase.verifyError(@() minnesotabvarm(2, 1, Y, Psi=[1 4], lambda2=0.5), ...
                "MATLAB:TooManyInputs");
            testCase.verifyError( ...
                @() minnesotabvarm(2, 1, Y, Type="inw", Psi=[1 4], lambda4=1), ...
                "MATLAB:TooManyInputs");
        end

        function unknownTypeMessageListsAliases(testCase)
            import matlab.unittest.constraints.ContainsSubstring

            message = minnesotaPriorClassesTest.unknownTypeMessage();

            testCase.verifyThat(message, ContainsSubstring("Available Type aliases:"));
            testCase.verifyThat(message, ...
                ContainsSubstring("mniw   (Matrix-Normal-Inverse-Wishart):"));
            testCase.verifyThat(message, ...
                ContainsSubstring("inw    (Independent Normal-Wishart):"));
            testCase.verifyThat(message, ...
                ContainsSubstring("normal (fixed-Sigma Normal):"));
            testCase.verifyThat(message, ...
                ContainsSubstring("Hyphens, underscores, and spaces are ignored"));
        end

        function rejectsMalformedType(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            testCase.verifyError(@() minnesotabvarm(2, 1, Y, Type=42), ...
                "minnesotabvarm:invalidType");
            testCase.verifyError(@() minnesotabvarm(2, 1, Y, Type="mniw", Type="inw"), ...
                "minnesotabvarm:repeatedType");
        end

        % ---- the shared base ---------------------------------------------

        function abstractBaseCannotBeInstantiated(testCase)
            testCase.verifyError(@() svar.minnesotabvarmBase(), "MATLAB:class:abstract");
        end

        function everyFamilySharesTheBase(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            for type = ["mniw" "inw" "normal"]
                prior = minnesotabvarm(2, 1, Y, Type=type, Psi=[1 4]);
                testCase.verifyTrue(isa(prior, "svar.minnesotabvarmBase"), ...
                    "Type=" + type + " should share the Minnesota base.");
            end
        end

        function everyFamilyStoresTheSampleAndPsi(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            for type = ["mniw" "inw" "normal"]
                prior = minnesotabvarm(2, 1, Y, Type=type, Psi=[1 4]);
                testCase.verifyEqual(prior.Y, Y, AbsTol=0);
                testCase.verifyEqual(prior.ResidualVariances, [1 4], AbsTol=0);
            end
        end

        function everyFamilyEstimatesPsiFromTheSample(testCase)
            Y = minnesotaPriorClassesTest.sampleData();
            expected = svar.estimateResidualVariances(Y, 1, Method="conditional");

            for type = ["mniw" "inw" "normal"]
                prior = minnesotabvarm(2, 1, Y, Type=type, Psi="conditional");
                testCase.verifyEqual(prior.ResidualVariances, expected, AbsTol=1e-12);
            end
        end

        function everyFamilyRejectsAForeignSample(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            for type = ["mniw" "inw" "normal"]
                prior = minnesotabvarm(2, 1, Y, Type=type, Psi=[1 4]);
                shortName = extractAfter(string(class(prior)), "svar.");
                testCase.verifyError(@() estimate(prior, Y(1:10,:)), ...
                    shortName + ":sampleMismatch");
            end
        end

        function everyFamilyValidatesTheSample(testCase)
            Y = minnesotaPriorClassesTest.sampleData();
            expected = ["minnesotamniwbvarm" "minnesotainwbvarm" "minnesotanbvarm"];
            types = ["mniw" "inw" "normal"];

            for k = 1:numel(types)
                testCase.verifyError( ...
                    @() minnesotabvarm(3, 1, Y, Type=types(k)), ...
                    expected(k) + ":invalidData", "Type=" + types(k));
            end
        end

        % ---- family-specific shapes --------------------------------------

        function inwCovarianceUsesFullCoefficientLayout(testCase)
            Y = minnesotaPriorClassesTest.sampleData(3);
            prior = svar.minnesotainwbvarm(3, 2, Y, Psi=[1 4 9]);
            numCoefficients = minnesotaPriorClassesTest.numEquationCoefficients(prior);

            testCase.verifySize(prior.Mu, [numCoefficients*prior.NumSeries 1]);
            testCase.verifySize(prior.V, ...
                [numCoefficients*prior.NumSeries numCoefficients*prior.NumSeries]);
            testCase.verifySize(prior.Omega, [prior.NumSeries prior.NumSeries]);
        end

        function normalCovarianceUsesFullCoefficientLayout(testCase)
            Y = minnesotaPriorClassesTest.sampleData(3);
            prior = svar.minnesotanbvarm(3, 2, Y, Psi=[1 4 9]);
            numCoefficients = minnesotaPriorClassesTest.numEquationCoefficients(prior);

            testCase.verifySize(prior.Mu, [numCoefficients*prior.NumSeries 1]);
            testCase.verifySize(prior.V, ...
                [numCoefficients*prior.NumSeries numCoefficients*prior.NumSeries]);
        end

        function normalPriorSetsFixedSigmaFromPsi(testCase)
            Y = minnesotaPriorClassesTest.sampleData(3);

            prior = svar.minnesotanbvarm(3, 2, Y, Psi=[1 4 9]);

            testCase.verifyEqual(prior.Sigma, diag([1 4 9]), AbsTol=0);
        end

        function normalPriorShrinksOnlyCrossLagCoefficientsWithLambda2(testCase)
            Y = minnesotaPriorClassesTest.sampleData();

            prior = svar.minnesotanbvarm(2, 1, Y, Psi=[1 4], ...
                lambda1=0.2, lambda2=0.5, lambda3=1);

            testCase.verifyEqual(diag(prior.V), ...
                [0.04; 0.0025; 1e4; 0.04; 0.04; 1e4], AbsTol=1e-14);
        end

        function inwAndNormalAgreeOnTheCoefficientCovariance(testCase)
            % Both families share buildIndependentCoefficientCovariance on the
            % base, so the same hyperparameters must give the same V.
            Y = minnesotaPriorClassesTest.sampleData();

            independent = svar.minnesotainwbvarm(2, 2, Y, Psi=[1 4], lambda2=0.3);
            normal = svar.minnesotanbvarm(2, 2, Y, Psi=[1 4], lambda2=0.3);

            testCase.verifyEqual(independent.V, normal.V, AbsTol=0);
        end

        function normalPriorEstimateReturnsNormalPosterior(testCase)
            Y = minnesotaPriorClassesTest.sampleData();
            prior = svar.minnesotanbvarm(2, 1, Y, Psi=[1 4], ...
                lambda1=0.2, lambda2=0.5, lambda3=1);

            testCase.verifyClass(estimate(prior, Display="off"), "normalbvarm");
        end

        function inwPriorEstimateReturnsEmpiricalPosterior(testCase)
            rng(0);
            Y = minnesotaPriorClassesTest.sampleData();
            prior = svar.minnesotainwbvarm(2, 1, Y, Psi=[1 4]);

            testCase.verifyClass(estimate(prior, Display="off"), "empiricalbvarm");
        end

        function normalPriorMarginalLikelihoodMethodsAgree(testCase)
            Y = minnesotaPriorClassesTest.sampleData();
            prior = svar.minnesotanbvarm(2, 1, Y, Psi=[1 4], ...
                lambda1=0.2, lambda2=0.5, lambda3=1);

            logML = logMarginalLikelihood(prior);
            ml = marginalLikelihood(prior);

            testCase.verifyTrue(isfinite(logML));
            testCase.verifyEqual(ml, exp(logML), AbsTol=1e-12);
        end

    end

    methods (Static, Access = private)
        function n = numEquationCoefficients(prior)
            n = prior.NumSeries*prior.P ...
                + double(prior.IncludeConstant) ...
                + double(prior.IncludeTrend) ...
                + prior.NumPredictors;
        end

        function message = unknownTypeMessage()
            try
                minnesotabvarm(2, 1, minnesotaPriorClassesTest.sampleData(), ...
                    Type="unknown");
                message = "";
            catch ME
                message = string(ME.message);
            end
        end

        function Y = sampleData(numseries)
            arguments
                numseries (1,1) double = 2
            end
            rng(0);
            Y = cumsum(randn(60, numseries))/10 + 1;
        end
    end
end
