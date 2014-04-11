function [pdfuncs,pdbranch,suc]=SetupMWPeriodDoubling(funcs,branch,ind,varargin)
%% initialize continuation of torus or period doubling bifurcations of periodic orbits
% simple wrapper to have a sensible name and set number of trivial Floquet
% multipliers to two
%
% $Id$
%
[pdfuncs,pdbranch,suc]=SetupMWTorusBifurcation(funcs,branch,ind,varargin{:});
end
