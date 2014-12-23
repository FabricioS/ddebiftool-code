function [hoho, success] = p_tohoho(point)

% function zeho = p_togenh(point)
% INPUT:
%	point: hopf point
% OUTPUT:
%	hoho: uncorrected starting guess for zero hopf point
%   success: whether conversion was successful

% (c) DDE-BIFTOOL v. 1.01, 14/07/2000

% Set success
success = 1;
hoho = point;

if strcmp(point.kind, 'hopf')
   if ~isfield(point,'stability') || ~isfield(point.stability,'l1') || isempty(point.stability.l1)
      error('P_TOHOHO: point does not contain stability information.');
   end
   hoho.kind = 'hoho';
   hoho.flag = '';
   if ~isfield(hoho,'nmfm')
      hoho.nmfm = [];
   end
   hoho.omega1 = point.omega;
   hoho = rmfield(hoho,'omega');
   roots = hoho.stability.l1;
   root1 = point.omega*sqrt(-1);
   root2 = -root1;
   roots = roots(abs(roots - root1) > 1e-8);
   roots = roots(abs(roots - root2) > 1e-8);
   imagpair = roots(abs(real(roots)) < 1e-6);
   if length(imagpair) < 2
      fprintf('P_TOHOHO: no second imaginary pair!');
      success = 0;
   end
   omega2 = abs(imag(imagpair(1)));
   hoho.omega2 = omega2;
else
   fprintf('P_TOHOHO: only hopf points can be converted into double hopf.\n');
   success = 0;
end

return