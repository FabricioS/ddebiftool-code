function p_splot(point)
%% plot spectrum of point with stability information
% function p_splot(point)
% INPUT:
%	point point whose stability needs plotting

% (c) DDE-BIFTOOL v. 1.02, 21/09/2001
%
% $Id$
%
%%
if strcmp(point.kind,'stst') || strcmp(point.kind,'fold') || strcmp(point.kind,'hopf')
  if isfield(point.stability,'l0')
    root_plt(point.stability.l0,point.stability.l1,point.stability.n1);
  end;
  xlabel('\Re(\lambda)');
  ylabel('\Im(\lambda)');
elseif strcmp(point.kind,'psol')
  if isfield(point.stability,'mu')
    mult_plt(point.stability.mu);
  end;
  xlabel('\Re(\mu)');
  ylabel('\Im(\mu)');
else
  err=point.kind;
  error('P_SPLOT: point kind %s not recognized.',err);
end;

end
