function [genh, success] = p_togenh(point)

% function genh = p_togenh(point)
% INPUT:
%	point: hopf point
% OUTPUT:
%	genh: uncorrected starting guess for generalized hopf point
%   success: whether conversion was successful

% (c) DDE-BIFTOOL v. 1.01, 14/07/2000

% Set success
success = 1;

genh = point;

if strcmp(point.kind, 'hopf')
   genh.kind = 'genh';
   genh.flag = '';
   if ~isfield(genh,'nmfm')
      genh.nmfm = [];
   end
   if ~isfield(genh.nmfm,'L2')
      genh.nmfm.L2 = NaN;
   end
else
   fprintf('P_TOGENH: only hopf points can be converted into generalized hopf.\n');
   success = 0;
end

return
