function fold=p_tofold(point)

% function fold_point=p_tofold(point)
% INPUT:
%	point point near fold point
% OUTPUT:
%	fold_point starting guess for fold point derived from point

% (c) DDE-BIFTOOL v. 1.00, 08/04/2000

fold.kind='fold';
fold.parameter=point.parameter;

if point.kind=='stst' | point.kind=='hopf'
  x=point.x;
elseif point.kind=='fold'
  error('P_TOFOLD: point is already fold.');
elseif point.kind=='psol'
  error('P_TOFOLD: periodic psol to fold not supported.');
end;

fold.x=x;
[D,dD]=root_cha(x,point.parameter,0);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2)));
fold.v=real(E1(:,i2));

return;
