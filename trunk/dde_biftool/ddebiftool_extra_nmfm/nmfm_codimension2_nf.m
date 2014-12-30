%% codimension-2 normal form for special point encountered along codimension-1 branch
%%
function [nf,nflow,br_ref,indbif]=...
    nmfm_codimension2_nf(funcs,branch,inds,detect,nmfm_compute,varargin)
%% Input
%
% * funcs: problem functions
% * branch: codim1 branch along which codim-2 poitn was encountered
% * inds: array of two successive indices bracing codimension-2 point
% * detect: monitoring function of type res=detect(p) for points of branch,
% or string ('hoho', 'zeho','hoze')
% * nmfm_compute: function of type newpoint=nmfm_type(funcs,point) where
% type is hoho, zeho etc
%
%% Output
%
% * nf: point with normal form
% * nflow: if numerical finite differences are used then computation is
% done twice, once with higher order, once with lower order, this output is
% the result with lower order. Use the difference between nf and nflow to
% estimate the error
% * br_ref: bisection refinements are performed along the branch before
% normal form calculation, br_ref is branch with additional points between
% indices inds
% * indbif: is the index in br_ref which is closest to the sign change of
% the monitoring function detect.
%
% $Id$
%
%%
default={'sys_mfderi',{}};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[br_ref,indbif]=br_bisection(funcs,branch,inds,detect,pass_on{:});
pt=br_ref.point(indbif);
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
if fin_diff
    newpoint=nmfm_compute(funcs1,pt);
    nf=newpoint;
    newpoint=nmfm_compute(funcs2,pt);
    nflow=newpoint;
else
    newpoint=nmfm_compute(funcs,pt);
    nf=newpoint;
    nflow=nf;
end
end
