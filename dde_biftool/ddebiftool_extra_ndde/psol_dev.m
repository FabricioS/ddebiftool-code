function [dpx]=psol_dev(profile1,t,x,m);

% function [dpx]=psol_dev(profile1,t,x,m);
% INPUT:
%	profile1 profile on mesh t
%	t representation points in [0,1]^(m*l+1)
%	x point(s) where to evaluate
%	m order polynomials
% OUTPUT: 
%	dpx value of profile derivative at x
% 
% (c) DDE-BIFTOOL v. 1.00, 08/04/2000
% Modified by David A.W. Barton 2006

for i=1:length(x)
    
  if (x(i) > 1.0) || (x(i) < 0.0)
      c=x(i)-floor(x(i));
  else
      c=x(i);
  end;

  index=length(t)-m;
  while (c<t(index)) 
    index=index-m; 
  end;

  Pa=poly_del(m,(c-t(index))/(t(index+m)-t(index)));  
  dpx(:,i)=profile1(:,index:index+m)*Pa'/(t(index+m)-t(index)); 

end;

return;