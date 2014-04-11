%% SetupPeriodDoubling - Initialize continuation of period doubling bifurcation
%%
function [pdfuncs,pdbranch,suc]=SetupPeriodDoubling(funcs,branch,ind,varargin)
%% 
% Simple wrapper around SetupTorusBifurcation to have a sensible name
% See <SetupTorusBifurcation.html> for description of input and output.
%
% <html>
% $Id$
% </html>
[pdfuncs,pdbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,varargin{:});
end
