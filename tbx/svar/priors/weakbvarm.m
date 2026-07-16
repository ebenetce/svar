classdef weakbvarm < conjugatebvarm
    %weakbvarm - Weak improper reduced-form Bayesian VAR prior
    %   PriorMdl = weakbvarm(NUMSERIES,NUMLAGS) creates a conjugate
    %   Bayesian VAR prior with diffuse coefficient uncertainty, zero
    %   inverse-Wishart scale, and zero prior degrees of freedom.
    %
    %   PriorMdl = weakbvarm(NUMSERIES,NUMLAGS,Name=VALUE) also
    %   specifies model options supported by conjugatebvarm, including
    %   SeriesNames, IncludeConstant, IncludeTrend, and NumPredictors.
    %
    %   This class implements the reduced-form weak prior used in Uhlig
    %   (2005), Appendix B: N0 = 0, n0 = 0, with S0 and B0 arbitrary.
    %   The posterior therefore has OLS coefficient mean, coefficient
    %   scale (X'*X)^(-1), inverse-Wishart scale equal to the residual
    %   sum of squares, and degrees of freedom equal to the effective
    %   sample size.
    %
    %   Posterior reduced-form parameter draws are available through the
    %   inherited simulate function, for example simulate(PriorMdl,Y)
    %   or simulate(PosteriorMdl).
    %
    %   weakbvarm functions:
    %       estimate  - Estimate the weak-prior reduced-form posterior
    %       simulate  - Simulate coefficients and innovations covariance matrix
    %       forecast  - Forecast responses from Bayesian vector autoregression (VAR) model
    %       simsmooth - Simulation smoother of Bayesian vector autoregression (VAR) model
    %       summarize - Distribution summary statistics of Bayesian vector autoregression (VAR) model
    %
    %   weakbvarm properties:
    %       Mu     - Diffuse coefficient prior mean
    %       V      - Diffuse coefficient prior scale
    %       Omega  - Zero inverse-Wishart scale matrix
    %       DoF    - Zero inverse-Wishart prior degrees of freedom
    %
    %   See also conjugatebvarm, diffusebvarm, uniformirbvarm

    methods
        function obj = weakbvarm(numseries, numlags, nvp)
            arguments
                numseries (1,1) double {mustBePositive, mustBeInteger}
                numlags (1,1) double {mustBeNonnegative, mustBeInteger}
                nvp.Description (1,1) string = string(missing)
                nvp.SeriesNames (1,:) string = string(missing)
                nvp.IncludeConstant (1,1) logical = true
                nvp.IncludeTrend (1,1) logical = false
                nvp.NumPredictors (1,1) double {mustBeNonnegative, mustBeInteger} = 0
            end

            % Series names
            if ismissing(nvp.SeriesNames)
                nvp = rmfield(nvp, "SeriesNames");
            end
            % Description
            if ismissing(nvp.Description)
                nvp = rmfield(nvp, "Description");
            end

            numDeterministics = nvp.IncludeConstant +  ...
                nvp.IncludeTrend + nvp.NumPredictors;
            k = numseries * numlags + numDeterministics;

            
            args = namedargs2cell(nvp);
            obj@conjugatebvarm(numseries, numlags, ...
                args{:}, ...
                Omega=zeros(numseries), ... % PphiBar
                DoF=0, ...                  % nnuBar
                Mu=zeros(k, numseries), ... % PpsiBar
                V=diag(inf(k, 1)));         % OomegaBar

            if ~isfield(nvp, "Description")
                obj.Description = obj.Description + " with weak improper NIW reduced-form prior";
            end
            
        end

        function [PosteriorMdl, Summary] = estimate(obj, varargin)
            [postBase, Summary] = estimate@conjugatebvarm(obj, varargin{:});

            % NB: avoids polymorphic estimate() which would return the
            % wrong class.
            PosteriorMdl = conjugatebvarm(postBase.NumSeries, postBase.P, ...
                SeriesNames=postBase.SeriesNames, ...
                IncludeConstant=postBase.IncludeConstant, ...
                IncludeTrend=postBase.IncludeTrend, ...
                NumPredictors=postBase.NumPredictors, ...
                Mu=postBase.Mu, ...
                V=postBase.V, ...
                Omega=postBase.Omega, ...
                DoF=postBase.DoF);
        end
    end
end
