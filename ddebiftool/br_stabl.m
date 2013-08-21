function branch=br_stabl(branch,skip,recompute)

% function st_branch=br_stabl(branch,skip,recompute)
% INPUT:
%	branch 
%	skip number of points to skip between stability computations
%	recompute if nonzero recompute stability info already present
% OUTPUT:
%	st_branch branch with stability information

% (c) DDE-BIFTOOL v. 1.02, 02/11/2000

ll=length(branch.point);

if ll<1 
  err=ll
  error('BR_STABL: branch is empty!');
end;

if ~isfield(branch.point(1),'stability')
  branch.point(1).stability=[];
end;

for i=1:skip+1:ll-1
  if isempty(branch.point(i).stability) | recompute
    branch.point(i).stability=p_stabil(branch.point(i),branch.method.stability);
  end;
end;

if isempty(branch.point(ll).stability) | recompute
  branch.point(ll).stability=p_stabil(branch.point(ll),branch.method.stability);
end;

return;