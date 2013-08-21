function p_splot(point)

% function p_splot(point)
% INPUT:
%	point point whose stability needs plotting

% (c) DDE-BIFTOOL v. 1.02, 21/09/2001

if point.kind=='stst' | point.kind=='fold' | point.kind=='hopf',
  if isfield(point.stability,'l0'),
    root_plt(point.stability.l0,point.stability.l1,point.stability.n1);
  end;
  xlabel('\Re(\lambda)');
  ylabel('\Im(\lambda)');
elseif point.kind=='psol',
  if isfield(point.stability,'mu'),
    mult_plt(point.stability.mu);
  end;
  xlabel('\Re(\mu)');
  ylabel('\Im(\mu)');
else
  err=point.kind
  error('P_SPLOT: point kind not recognized.');
end;

return;
