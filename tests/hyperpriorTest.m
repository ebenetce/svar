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

        function nativeParameterizationStoresParamsUnconverted(testCase)
            % The GLP (2012) hyperprior on the residual-variance scale:
            % Inverse-Gamma with shape = scale = 0.02^2, which has neither a
            % mean nor a variance and so cannot be reached from (mode, sd).
            hp = hyperprior("InverseGamma", 0.02^2, 0.02^2, ...
                Parameterization="native", Bounds=[1e-5 1], X0=1e-3);

            testCase.verifyEqual(hp.Params, [0.02^2 0.02^2], AbsTol=0);
            testCase.verifyEqual(hp.meanOf, Inf);
            testCase.verifyTrue(isfinite(hp.logpdf(hp.X0)));
        end

        function nativeLogpdfMatchesInverseGammaDensity(testCase)
            a = 0.02^2;
            hp = hyperprior("InverseGamma", a, a, Parameterization="native", ...
                Bounds=[1e-5 1], X0=1e-3);
            x = [1e-4 1e-3 1e-2];

            expected = a*log(a) - gammaln(a) - (a + 1)*log(x) - a./x;
            testCase.verifyEqual(hp.logpdf(x), expected, RelTol=1e-12);
        end

        function momentParameterizationIsUnchangedByDefault(testCase)
            % Same call with and without the explicit default must agree.
            implicitDefault = hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5]);
            explicitDefault = hyperprior("Gamma", 0.2, 0.4, Bounds=[1e-4 5], ...
                Parameterization="moments");

            testCase.verifyEqual(explicitDefault.Params, ...
                implicitDefault.Params, AbsTol=0);
        end

        function diffuseNativePriorRejectsDefaultBounds(testCase)
            % The quantile box is unusable here, and must be reported rather
            % than silently stored - GAMINV does not converge that far out.
            testCase.verifyError( ...
                @() hyperprior("InverseGamma", 0.02^2, 0.02^2, ...
                    Parameterization="native"), ...
                "hyperprior:unusableDefaultBounds");
        end

        function unknownParameterizationIsRejected(testCase)
            testCase.verifyError( ...
                @() hyperprior("Gamma", 0.2, 0.4, Parameterization="native2"), ...
                "MATLAB:validators:mustBeMember");
        end
    end
end
