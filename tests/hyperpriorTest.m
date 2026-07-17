classdef hyperpriorTest < matlab.unittest.TestCase
    %hyperpriorTest Tests for hyperprior value objects.

    methods (Test)
        function explicitBoundsAndX0AreStored(testCase)
            hp = hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5], X0=0.3);

            testCase.verifyEqual(hp.Bounds, [1e-4 5], AbsTol=0);
            testCase.verifyEqual(hp.X0, 0.3, AbsTol=0);
        end

        function defaultBoundsAndX0AreFinite(testCase)
            hp = hyperprior("Gamma", 0.2, 0.4);

            testCase.verifySize(hp.Bounds, [1 2]);
            testCase.verifyTrue(all(isfinite(hp.Bounds)));
            testCase.verifyTrue(isfinite(hp.X0));
        end

        function logpdfIsFiniteAtX0(testCase)
            hp = hyperprior("InverseGamma", 0.5, 0.25, Bounds=[0.05 5], X0=0.5);

            testCase.verifyTrue(isfinite(hp.logpdf(hp.X0)));
        end

        function invalidBoundsError(testCase)
            testCase.verifyError( ...
                @() hyperprior("Gamma", 0.2, 0.4, Bounds=[5 1]), ...
                "hyperprior:invalidBounds");
        end

        function positiveSupportDistributionsAcceptInfiniteUpperBound(testCase)
            gammaPrior = hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 Inf]);
            inverseGammaPrior = hyperprior("InverseGamma", 0.2, 0.4, ...
                Bounds=[1e-4 Inf]);

            testCase.verifyEqual(gammaPrior.Bounds, [1e-4 Inf], AbsTol=0);
            testCase.verifyEqual(inverseGammaPrior.Bounds, [1e-4 Inf], AbsTol=0);
        end

        function betaRejectsBoundsOutsideSupport(testCase)
            testCase.verifyError( ...
                @() hyperprior("Beta", 0.5, 0.1, Bounds=[1e-4 Inf]), ...
                "hyperprior:invalidBounds");
        end

        function invalidX0Error(testCase)
            testCase.verifyError( ...
                @() hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5], X0=10), ...
                "hyperprior:invalidX0");
        end
    end
end
