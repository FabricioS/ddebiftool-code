function [t_new]=psol_msh(t,m,profile,l_t_new,m_new);

% function [t_new]=psol_msh(t,m,profile,l_t_new,m_new);
% INPUT:
%	t representation points
%	m number of collocation points (per interval)
%	profile solution profile
%	l_t_new size of new mesh (number of intervals)
%	m_new new number of collocation points (per interval)
% OUTPUT:
%	t_new new, adapted mesh

% (c) DDE-BIFTOOL v. 1.00, 11/03/2000

n=size(profile,1);
l=(length(t)-1)/m;

if l~=floor(l),
  err=[length(t) m]
  error('PSOL_MSH: t does not contain l intervals of m points!');
end;

ti=t(1:m:length(t));
dti=ti(2:length(ti))-ti(1:length(ti)-1);

[ti_new]=auto_msh(n,profile,l,m,ti,dti,l_t_new);

for i=1:l_t_new,
  t_new(m_new*(i-1)+1:m_new*i)=ti_new(i)+(ti_new(i+1)-ti_new(i))*(0:m_new-1)/m_new;
end;

t_new(length(t_new)+1)=1;

return;