classdef historicalDecompositionTest < matlab.unittest.TestCase
    %historicalDecompositionTest Tests for historicalDecomposition input handling.

    methods (Test)
        function acceptsVarmDirectly(testCase)
            % A fully specified varm is the natural input and must not need
            % conversion (regression: the function used to force bvar2var).
            [mdl, Y] = historicalDecompositionTest.stableModelAndData();
            impact = chol(mdl.Covariance, "lower");

            HD = historicalDecomposition(mdl, impact, Y);

            testCase.verifySize(HD.Contributions, ...
                [size(Y,1) - mdl.P, mdl.NumSeries, mdl.NumSeries]);
            testCase.verifyEqual(HD.Model, mdl);
        end

        function varmAndBvarInputsAgree(testCase)
            % Passing a BVAR object and its bvar2var conversion must give an
            % identical decomposition - the whole point of accepting both.
            [~, Y] = historicalDecompositionTest.stableModelAndData();
            prior = minnesotabvarm(size(Y,2), 2, Y);
            post  = estimate(prior, Display="off");
            varMdl = bvar2var(post);
            impact = chol(varMdl.Covariance, "lower");

            fromBvar = historicalDecomposition(post,   impact, Y);
            fromVarm = historicalDecomposition(varMdl, impact, Y);

            testCase.verifyEqual(fromVarm.Contributions, fromBvar.Contributions, ...
                AbsTol=1e-12);
            testCase.verifyEqual(fromVarm.StructuralShocks, fromBvar.StructuralShocks, ...
                AbsTol=1e-12);
        end

        function rejectsNonModelInput(testCase)
            [~, Y] = historicalDecompositionTest.stableModelAndData();

            testCase.verifyError( ...
                @() historicalDecomposition(struct("NumSeries", 2), eye(2), Y), ...
                "historicalDecomposition:InvalidModel");
        end

        function validatesImpactSize(testCase)
            [mdl, Y] = historicalDecompositionTest.stableModelAndData();

            testCase.verifyError( ...
                @() historicalDecomposition(mdl, eye(mdl.NumSeries + 1), Y), ...
                "historicalDecomposition:InvalidImpactSize");
        end
    end

    methods (Static, Access = private)
        function [mdl, Y] = stableModelAndData()
            mdl = varm(2, 2);
            mdl.Constant = [0.1; -0.1];
            mdl.AR = {[0.5 0.1; 0.0 0.4], [0.1 0.0; 0.1 0.2]};
            mdl.Covariance = [1.0 0.3; 0.3 0.8];
            rng(42, "twister");
            Y = simulate(mdl, 120);
        end
    end
end
