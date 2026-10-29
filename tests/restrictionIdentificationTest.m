classdef restrictionIdentificationTest < matlab.unittest.TestCase

    methods (TestMethodSetup)
        function seedRandomStream(testCase)
            originalStream = rng;
            testCase.addTeardown(@() rng(originalStream));
            rng(7, "twister");
        end
    end

    methods (Test)
        function storesNamedSignRestriction(testCase)
            restriction = svar.SignRestriction(1, ...
                Response="Output", Shock="Demand", Horizon=0:2);

            testCase.verifyEqual(restriction.Response, "Output");
            testCase.verifyEqual(restriction.Shock, "Demand");
            testCase.verifyEqual(restriction.requiredHorizon(), 2);
        end

        function storesNamedZeroRestriction(testCase)
            restriction = svar.ZeroRestriction( ...
                Response="Output", Shock="Demand");

            testCase.verifyEqual(restriction.Response, "Output");
            testCase.verifyEqual(restriction.Shock, "Demand");
            testCase.verifyEqual(restriction.Horizon, 0);
        end

        function drawsAdmissibleRotationsForFittedVar(testCase)
            mdl = restrictionIdentificationTest.model();
            restrictions = [ ...
                svar.SignRestriction(1, ...
                    Response="Output", Shock="Demand", Horizon=0), ...
                svar.SignRestriction(1, ...
                    Response="Prices", Shock="Demand", Horizon=0)];

            [impact, info] = svar.identify(mdl, restrictions, ...
                ShockNames=["Demand" "Supply"], NumDraws=4);

            testCase.verifySize(impact, [2 2 4]);
            testCase.verifyEqual(info.Impact, impact, AbsTol=0);
            testCase.verifyEqual(info.Method, "independent");
            testCase.verifyEqual(info.NumDraws, 4);
            testCase.verifyGreaterThan(info.NumCandidates, 3);
            testCase.verifyGreaterThanOrEqual(info.AcceptanceRate, 0);
            testCase.verifyGreaterThanOrEqual(squeeze(impact(1, 1, :)), 0);
            testCase.verifyGreaterThanOrEqual(squeeze(impact(2, 1, :)), 0);
        end

        function constructsImpactZeroRestrictions(testCase)
            mdl = restrictionIdentificationTest.model();
            restrictions = [ ...
                svar.ZeroRestriction(Response="Output", Shock="Demand"), ...
                svar.SignRestriction(1, ...
                    Response="Prices", Shock="Demand", Horizon=0)];

            [impact, info] = svar.identify(mdl, restrictions, ...
                ShockNames=["Demand" "Supply"], NumDraws=4);

            testCase.verifyEqual(squeeze(impact(1, 1, :)), zeros(4, 1), ...
                AbsTol=1e-12);
            testCase.verifyGreaterThanOrEqual(squeeze(impact(2, 1, :)), 0);
            testCase.verifyEqual( ...
                pagemtimes(impact, "none", permute(impact, [2 1 3]), "none"), ...
                repmat(mdl.Covariance, 1, 1, 4), AbsTol=1e-12);
            testCase.verifyTrue(info.Restrictions.hasZeroRestrictions());
        end

        function constructsZeroForLaterShock(testCase)
            mdl = restrictionIdentificationTest.model();
            restrictions = [ ...
                svar.ZeroRestriction(Response="Output", Shock="Supply"), ...
                svar.SignRestriction(1, ...
                    Response="Prices", Shock="Supply", Horizon=0)];

            impact = svar.identify(mdl, restrictions, ...
                ShockNames=["Demand" "Supply"], NumDraws=4);

            testCase.verifyEqual(squeeze(impact(1, 2, :)), zeros(4, 1), ...
                AbsTol=1e-12);
            testCase.verifyGreaterThanOrEqual(squeeze(impact(2, 2, :)), 0);
        end

        function rejectsNonimpactZeroRestriction(testCase)
            testCase.verifyError(@() svar.ZeroRestriction( ...
                Response="Output", Shock="Demand", Horizon=1), ...
                "svar:ZeroRestriction:UnsupportedHorizon");
        end

        function retainsPosteriorDrawMetadata(testCase)
            Y = randn(50, 2);
            prior = minnesotabvarm(2, 1, Y, ...
                SeriesNames=["Output" "Prices"]);
            posterior = estimate(prior, Display="off");
            restrictions = svar.SignRestriction(1, ...
                Response="Output", Shock="Demand", Horizon=0:1);

            [impact, info] = svar.identify(posterior, restrictions, ...
                ShockNames=["Demand" "Supply"], NumDraws=2);

            testCase.verifySize(impact, [2 2 2]);
            testCase.verifySize(info.Coefficients, [numel(posterior.Mu) 2]);
            testCase.verifySize(info.Sigma, [2 2 2]);
            testCase.verifySize(info.Rotation, [2 2 2]);
            testCase.verifyEqual( ...
                pagemtimes(impact, "none", permute(impact, [2 1 3]), "none"), ...
                info.Sigma, AbsTol=1e-10);
        end

        function drawsSignRestrictedEllipticalChain(testCase)
            Y = randn(50, 2);
            prior = minnesotabvarm(2, 1, Y, ...
                SeriesNames=["Output" "Prices"]);
            posterior = estimate(prior, Display="off");
            restrictions = svar.SignRestriction(1, ...
                Response="Output", Shock="Demand", Horizon=0);

            [impact, info] = svar.identify(posterior, restrictions, ...
                ShockNames=["Demand" "Supply"], Sampler="elliptical", ...
                NumDraws=3, BurnIn=2, NumParticles=20);

            testCase.verifyEqual(info.Method, "elliptical");
            testCase.verifyEqual(info.InitializationMethod, "smc");
            testCase.verifySize(impact, [2 2 3]);
            testCase.verifyGreaterThanOrEqual(squeeze(impact(1, 1, :)), 0);
            testCase.verifyEqual( ...
                pagemtimes(impact, "none", permute(impact, [2 1 3]), "none"), ...
                info.Sigma, AbsTol=1e-10);
        end

        function rejectsZeroRestrictionsForEllipticalSampler(testCase)
            Y = randn(50, 2);
            prior = minnesotabvarm(2, 1, Y, ...
                SeriesNames=["Output" "Prices"]);
            posterior = estimate(prior, Display="off");
            restrictions = svar.ZeroRestriction( ...
                Response="Output", Shock="Demand");

            testCase.verifyError(@() svar.identify(posterior, restrictions, ...
                ShockNames=["Demand" "Supply"], Sampler="elliptical"), ...
                "svar:identify:EllipticalZeroRestrictions");
        end

        function rejectsUnknownNamedResponse(testCase)
            mdl = restrictionIdentificationTest.model();
            restrictions = svar.SignRestriction(1, ...
                Response="Missing", Shock="Demand", Horizon=0);

            testCase.verifyError(@() svar.identify(mdl, restrictions, ...
                ShockNames=["Demand" "Supply"]), ...
                "svar:SignRestriction:UnknownResponse");
        end

        function rejectsIncompleteShockNames(testCase)
            mdl = restrictionIdentificationTest.model();
            restrictions = svar.SignRestriction(1, ...
                Response="Output", Shock="Demand", Horizon=0);

            testCase.verifyError(@() svar.identify(mdl, restrictions, ...
                ShockNames="Demand"), "svar:identify:InvalidShockNames");
        end
    end

    methods (Static, Access = private)
        function mdl = model()
            mdl = varm(2, 1);
            mdl.SeriesNames = ["Output" "Prices"];
            mdl.AR = {0.1*eye(2)};
            mdl.Covariance = [1 0.1; 0.1 1];
        end
    end
end
