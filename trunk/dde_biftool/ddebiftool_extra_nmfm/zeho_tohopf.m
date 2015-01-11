function hopf=zeho_tohopf(funcs,point,freqs) %#ok<INUSL,INUSD>
%% convert zero-Hopf point to Hopf bifurcation point
% function hopf_point=zeho_tohopf(funcs,point {,freqs})
% INPUT:
%   funcs problem functions
%	point with stability information 
%   optional freqs: frequency to be excluded from consideration
% OUTPUT:
%	hopf_point: hopf point
%
% $Id$
%
%%
hopf = struct(...
    'kind','hopf',...
    'parameter',point.parameter,...
    'x',point.x,...
    'v',point.v,...
    'omega',point.omega);
if isfield(point,'stability');
    hopf.stability=point.stability;
end
end