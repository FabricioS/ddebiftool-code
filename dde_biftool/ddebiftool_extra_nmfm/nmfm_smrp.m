function [smallest_real_part, stability,selectedroot] = nmfm_smrp(...
    funcs,point, stmethod,rmomega,threshold)
%% Compute smallest real part of imaginary pairs or of real eigenvalues
%
% * rmomega (logical): remove known imaginary roots
% * threshold (function): threshold(abs(imag(roots)) selects roots to
% consider
%
% $Id$
%
%%
if ~isfield(point, 'stability') || isempty(point.stability) || isempty(point.stability.l0)
	point.stability = p_stabil(funcs,point, stmethod);
end

stability = point.stability;
roots = stability.l0;

%% Remove roots closest to known eigenvalue pair
if rmomega
   [dum,ix]=sort(abs(roots - 1i*point.omega)); %#ok<ASGLU>
   roots(ix(1))=[];
   [dum,ix]=sort(abs(roots + 1i*point.omega)); %#ok<ASGLU>
   roots(ix(1))=[];
end

selectedroots = roots(threshold(abs(imag(roots))));
if isempty(selectedroots)
   smallest_real_part = 1e6;
   return
end
%% Assume all imaginary eigenvalues come in pairs
realparts = real(selectedroots);
[dum, rpind] = sort(abs(realparts)); %#ok<ASGLU>
smallest_real_part = realparts(rpind(1));
selectedroot=selectedroots(rpind(1));
end
