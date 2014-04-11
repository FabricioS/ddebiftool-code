function [trfuncs,trbranch,suc]=SetupMWTorusBifurcation(funcs,branch,ind,varargin)
%% initialize continuation of torus or period doubling bifurcations of periodic orbits
% simple wrapper to have a sensible name and set number of trivial Floquet
% multipliers to two
%
% $Id$
%
[trfuncs,trbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,'nremove',2,varargin{:});
end
