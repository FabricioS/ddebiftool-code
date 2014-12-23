function flaggedbranch = br_flag(branch)

% function flaggedbranch = br_flag(branch)
% Purpose:
%   Sets flag = '' on all points of the branch.
% INPUT:
%	branch 
% OUTPUT:
%	flaggedbranch

% (c) DDE-BIFTOOL v. 2.00, 23/12/2000

flaggedbranch = branch;
ll=length(branch.point);

if ll<1
   error('BR_FLAG: branch is empty!');
end;

for i=1:ll
   if ~isfield(flaggedbranch.point(i),'flag') || isempty(flaggedbranch.point(i).flag)
      flaggedbranch.point(i).flag = '';
   end
   if ~isfield(flaggedbranch.point(i),'nmfm')
      flaggedbranch.point(i).nmfm = [];
   end
   if ~isfield(flaggedbranch.point(i),'nvec')
      flaggedbranch.point(i).nvec = [];
   end
end;



return;
