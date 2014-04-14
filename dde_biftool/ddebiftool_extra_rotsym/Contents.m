% DDEBIFTOOL_EXTRA_ROTSYM
%
% Files
%   rot_cond                - extra phase condition needed to fix rotation phase
%   rot_deriv               - adapt user provided derivative to rotating coordinates
%   rot_rhs                 - right-hand side in rotating coordinates
%   set_rotfuncs            - fill in funcs structure for use with DDE-Biftool, rotational symmetric case
%   SetupMWFold             - initialize continuation of folds of modulated waves
%   SetupMWPeriodDoubling   - initialize continuation of torus or period doubling bifurcations of periodic orbits
%   SetupMWTorusBifurcation - initialize continuation of torus or period doubling bifurcations of periodic orbits
%   SetupRWFold             - initialize continuation of folds of periodic orbits
%   SetupRWHopf             - Initialize continuation of Hopf bifurcations of relative equilibria
%   sys_cond_MWFold         - constraints used for extended DDE in fold continuation for relative periodic orbits
%   sys_cond_RWFold         - constraints used for extended DDE in fold continuation of relative equilibria
%   sys_rhs_MWFold          - rhs of extended DDE for fold of modulated waves
%   sys_rhs_RWFold          - rhs of extended DDE for fold of relative equilibria
