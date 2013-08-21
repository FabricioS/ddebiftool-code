function p=poly_lgr(t,c);

% function p=poly_lgr(t,c);
% INPUT:
%	t lagrange points in R^m+1
% 	c evaluation point in R 
% OUTPUT:
%	p values of lagrange polynomials through t at c

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

m=length(t)-1;

% compute p:

for j=1:m+1
  p(j)=1;
  for k=1:m+1
    if k~=j
      p(j)=p(j)*(c-t(k))/(t(j)-t(k));
    end;
  end;
end;

return;
