%% Initialize continuation of Fold bifurcations
%%
function [foldbranch,suc]=SetupFold(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of stst along which Fold was discovered
% * |ind|: index of approximate Fold point
%
% Important optional inputs (name-value pairs)
% outputs
% foldbranch: Fold branch with first point (or two points)
% suc: flag whether corection was successful
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of additional continuation parameters
%  (in addition to free pars of branch
% * |'correc'| (logical, default true): apply |p_correc| to first points on fold
%   branch
% * |'dir'| (integer, default |[]|): which parameter to vary initially along fold
%   branch (foldbranch has only single point if dir is empty)
% * |'step'| (real, default |1e-3|): size of initial step if dir is non-empty
%
% All other named arguments are passed on to fields of |foldbranch|
%% Outputs
% 
% * |foldbranch|: Hopf branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
% $Id$
%
%% process options
default={'contpar',[],'correc',true,'dir',[],'step',1e-3};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% initialize branch of folds (pbranch)
foldbranch=df_brnch(funcs,options.contpar,'fold');
foldbranch.method.point=branch.method.point;
foldbranch.method.continuation=dde_set_options(foldbranch.method.continuation,pass_on,'pass_on');
foldbranch.method.point=dde_set_options(foldbranch.method.point,pass_on,'pass_on');
foldbranch.method.stability=dde_set_options(foldbranch.method.stability,pass_on,'pass_on');
if isempty(options.contpar)
    % if no continuation parameters are given use free parameters of branch
    foldbranch.parameter.free=branch.parameter.free;
elseif length(options.contpar)==1;
    % if single continuation parameter is given prepend to free parameters of branch
    foldbranch.parameter.free=[options.contpar,branch.parameter.free];
else
    foldbranch.parameter.free=options.contpar(:)';
end
foldbranch.parameter=dde_set_options(foldbranch.parameter,pass_on,'pass_on');
point=branch.point(ind);
if ~isfield(point,'stability')
    point.stability=p_stabil(funcs,point,branch.method.stability);
end
%% set up functions of extended system
%% create initial guess for correction
foldini0=p_tofold(funcs,point);
%% correct and add 2nd point if desired
[foldbranch,suc]=correct_ini(funcs,foldbranch,foldini0,...
    options.dir,options.step,options.correc);
end
