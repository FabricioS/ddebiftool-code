function [hbranch,suc]=SetupRWHopf(funcs,branch,ind,varargin)
%% Initialize continuation of Hopf bifurcations of relative equilibria
% simple wrapper to have a sensible name and set number of trivial Floquet
% multipliers to two
%
% $Id$
%
[hbranch,suc]=SetupHopf(funcs,branch,ind,varargin{:});
end
