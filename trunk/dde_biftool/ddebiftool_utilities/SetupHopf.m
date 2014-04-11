%% Initialize continuation of Hopf bifurcations
%%
function [hbranch,suc]=SetupHopf(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of stst along which Hopf was discovered
% * |ind|: index of approximate Hopf point
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of additional continuation parameters
%  (in addition to free pars of branch
% * |'correc'| (logical, default true): apply |p_correc| to first points on
% hopf branch
% * |'dir'| (integer, default |[]|): which parameter to vary initially
% along Hopf branch (hbranch has only single point if dir is empty)
% * |'step'| (real, default |1e-3|): size of initial step if dir is non-empty
%
% All other named arguments are passed on to fields of |hbranch|
%% Outputs
% 
% * |hbranch|: Fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%
% $Id$
%
%% process options
default={'contpar',[],'correc',true,'dir',[],'step',1e-3,'excludefreqs',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% initialize branch of folds (pbranch)
hbranch=branch;
hbranch=replace_branch_pars(hbranch,options.contpar,pass_on);
point=branch.point(ind);
if ~isfield(point,'stability') || isempty(point.stability)
    point.stability=p_stabil(funcs,point,branch.method.stability);
end
%% create initial guess for correction
hini0=p_tohopf(funcs,point,options.excludefreqs);
%% correct and add 2nd point if desired
[hbranch,suc]=correct_ini(funcs,hbranch,hini0,...
    options.dir,options.step,options.correc);
end
