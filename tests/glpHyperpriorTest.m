classdef glpHyperpriorTest < matlab.unittest.TestCase
    %glpHyperpriorTest Tests for GLP hyperprior inputs.

    methods (Test)
        function acceptsFixedNumericPsi(testCase)
            Y = glpHyperpriorTest.sampleData();
            opts = glpHyperpriorTest.optimOptions();

            [mdl, info] = glp(2, 1, Y, [0.5 1.0], ...
                OptimOptions=opts);

            testCase.verifyClass(mdl, "minnesotamniwbvarm");
            testCase.verifyFalse(info.PsiFree);
            testCase.verifyEqual(info.FinalPsi, [0.5 1.0], AbsTol=0);
        end

        function acceptsNumericPsiBounds(testCase)
            Y = glpHyperpriorTest.sampleData();
            opts = glpHyperpriorTest.optimOptions();

            [mdl, info] = glp(2, 1, Y, [0.1 0.2; 2.0 3.0], ...
                OptimOptions=opts);

            testCase.verifyClass(mdl, "minnesotamniwbvarm");
            testCase.verifyTrue(info.PsiFree);
            testCase.verifyEqual(info.PsiNames, ["Psi1" "Psi2"]);
        end

        function acceptsPsiHyperpriorArray(testCase)
            Y = glpHyperpriorTest.sampleData();
            opts = glpHyperpriorTest.optimOptions();
            spec = minnesotaSpec("mniw", ...
                lambda1=hyperprior("Gamma", 0.2, 0.4, Bounds=[0.05 0.5]));
            Psi = arrayfun(@(x) hyperprior("InverseGamma", x, 0.5*x, ...
                X0=x, Bounds=[0.1 10]*x), [0.5 1.0]);

            [mdl, info] = glp(2, 1, Y, Psi, ...
                Spec=spec, OptimOptions=opts);

            testCase.verifyClass(mdl, "minnesotamniwbvarm");
            testCase.verifyTrue(info.PsiFree);
            testCase.verifyTrue(info.UsedHyperprior);
            testCase.verifyEqual(numel(info.X0), 3);
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

        function opts = optimOptions()
            opts = optimoptions("fmincon", ...
                Display="off", ...
                MaxIterations=1, ...
                MaxFunctionEvaluations=30);
        end
    end
end
