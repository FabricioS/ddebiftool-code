%% first Lyapunov coefficient L1 in Hopf point along stst branch or Hopf branch
function [L1,L1low]=HopfLyapunovCoefficients(funcs,branch,varargin)
%% Input
%
% * funcs: problem functions
% * branch: hopf branch or array of points of kind 'hopf'
% * optional name-value pairs: 
% * 'use_nullvectors' (boolean, default false) use nullvectors from
% previous point along branch to compute new nullvectors by bordering
% (otherwise, svd is used)
% * 'sys_mfderi': update sys_mfderi field of funcs (useful for
% finite-difference approzimation only)
%
%% Output
%
% * L1: array of L1 coefficients along branch
% * L1low: equals L1 if funcs.sys_mfderi_provided is true, otherwise (for
% finite-difference approximation) a low-er-order estimate of L1
% coefficients. Use L1-L1low to estimate the accuracy of L1. If this is
% large, finite-differences are problematic. If this is small it may(!) be
% ok.
%
% $Id$
%
%%
default={'use_nullvectors',false,'sys_mfderi',{}};
options=dde_set_options(default,varargin);
if isfield(branch,'point')
    pt=branch.point;
else
    pt=branch;
end
npt=length(pt);
if funcs.sys_mfderi_provided
    fin_diff=false;
else
    fin_diff=true;
    if ~isempty(options.sys_mfderi)
        funcs1=set_funcs(funcs,'sys_mfderi',options.sys_mfderi);
    else
        funcs1=funcs;
    end
    funcs2=set_funcs(funcs1,'sys_mfderi',[options.sys_mfderi,{'output',2}]);
end
nullpoint={};
L1=NaN(1,npt);
L1low=L1;
for i=1:npt
    currpt=pt(i);
    if fin_diff
        newpoint=nmfm_hopf(funcs1,currpt,nullpoint{:});
        L1(i)=newpoint.nmfm.L1;
        newpoint=nmfm_hopf(funcs2,currpt,nullpoint{:});
        L1low(i)=newpoint.nmfm.L1;
    else
        newpoint=nmfm_hopf(funcs, currpt,nullpoint{:});
        L1(i)=newpoint.nmfm.L1;
        L1low(i)=L1(i);
    end
    if options.use_nullvectors
        nullpoint={newpoint};
    end
end
end
