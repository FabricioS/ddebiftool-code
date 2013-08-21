function [px]=psol_eva(profile1,t,x,m);

% function [px]=psol_eva(profile1,t,x,m);
% INPUT:
%	profile1 profile on mesh t
%	t representation points in [0,1]^(m*l+1)
%	x point(s) where to evaluate
%	m order polynomials
% OUTPUT: 
%	px value of profile at x

% (c) DDE-BIFTOOL v. 1.00, 08/04/2000

for i=1:length(x)

  c=x(i)-floor(x(i));

  index=length(t)-m;
  while (c<t(index)) 
    index=index-m; 
  end;

  Pa=poly_elg(m,(c-t(index))/(t(index+m)-t(index)));  
  px(:,i)=profile1(:,index:index+m)*Pa'; 

end;

return;