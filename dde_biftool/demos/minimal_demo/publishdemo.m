%%  DDE-BIFTOOL minimal_demo - Duffing oscillator with delayed feedback
% run all scripts of this for testing
% This script runs all demos files in a single run for testing purposes.
% See <html/minimal_demo.html> for the published version of the single scripts
%
% <html>
% $Id$
% </html>
%% Description, load path and define system
publish('minimal_demo');
%% Steady-state bifurcations and periodic orbits
opts={'maxOutputLines',20};
publish('minimal_demo_stst_psol',opts{:});
%% local bifurcations of periodic orbits
publish('minimal_demo_extra_psol',opts{:});
