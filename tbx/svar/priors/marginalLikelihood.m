function ml = marginalLikelihood(Mdl, Y, opts)
%MARGINALLIKELIHOOD Marginal likelihood for supported BVAR priors.
%   ML = MARGINALLIKELIHOOD(MDL,Y) returns p(Y | MDL) for supported BVAR
%   prior objects by exponentiating LOGMARGINALLIKELIHOOD(MDL,Y).
%
%   ML = MARGINALLIKELIHOOD(MDL,Y,X=X,Y0=Y0) supplies exogenous predictors
%   or presample responses to the sufficient-statistics calculation.
%
%   [ML,DETAILS] = MARGINALLIKELIHOOD(___) also returns the details from the
%   log marginal likelihood calculation.
%
%   See also LOGMARGINALLIKELIHOOD

arguments
    Mdl (1,1)
    Y {mustBeNonempty}
    opts.X = []
    opts.Y0 = []
end

logML = logMarginalLikelihood(Mdl,  Y,  X=opts.X, Y0=opts.Y0);
ml = exp(logML);

end
