function [smallest_real_part, fullroots] = nmfm_smrpip(funcs,point, stmethod, rmomega)
%% Compute smallest real part of imaginary pairs
%
% $Id$
%
%%
if ~isfield(point, 'stability') || isempty(point.stability) || isempty(point.stability.l1)
	point.stability = p_stabil(funcs,point, stmethod);
end

roots = point.stability.l1;
fullroots = roots;

% Remove known eigenvalue pair
if rmomega
   root1 = point.omega*sqrt(-1);
   root2 = -root1;
   roots = roots(abs(roots - root1) > 1e-8);
   roots = roots(abs(roots - root2) > 1e-8);
end

imagroots = roots(abs(imag(roots)) > 1e-6);
if isempty(imagroots)
   smallest_real_part = 1e6;
   return;
end
% Assume all imaginary eigenvalues come in pairs
realparts = real(imagroots);
[~, rpind] = sort(abs(realparts));
realparts = realparts(rpind);
smallest_real_part = realparts(1);

end
