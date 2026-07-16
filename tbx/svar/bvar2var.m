function varMdl = bvar2var(bvarMdl)
%bvar2var  Convert a bvar model to a varm model
%   varMdl = bvar2var(bvarMdl) Converts a Bayesian VAR model object into an
%   equivalent varm model. The function maps the bvar model properties
%   (Constant, AR, Trend, Beta, Covariance, SeriesNames) into a varm object
%   and returns that varm model as the output.
%
%   Examples:
%       bvarMdl = conjugatebvarm(3, 2); 
%       varMdl = bvar2var(bvarMdl);
%
%       % If a bvar model already exists in workspace
%       load('myBvarModel.mat','bvarMdl'); 
%       varMdl = bvar2var(bvarMdl);
%
%   See also VARM, BVAR

arguments
    bvarMdl (1,1) {mustBeVarOrBvar(bvarMdl)}
end

if isa(bvarMdl, 'varm')
    varMdl = bvarMdl;
    return
end

varMdl = varm(AR=bvarMdl.AR, SeriesNames=bvarMdl.SeriesNames, ...
    Covariance=bvarMdl.Covariance);

if ~isempty(bvarMdl.Constant)
    varMdl.Constant = bvarMdl.Constant;
end

if ~isempty(bvarMdl.Trend)
    varMdl.Trend = bvarMdl.Trend;
end

if ~isempty(bvarMdl.Beta)
    varMdl.Beta = bvarMdl.Beta;
end

end
