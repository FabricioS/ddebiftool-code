function dp=poly_dla(t,c);

% function dp=poly_dla(t,c);
% INPUT:
%       t lagrange points in R^m+1
%       c evaluation point in R 
% OUTPUT:
%       dp values of derivative of lagrange polynomials through t at c

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

m=length(t)-1;

% compute dp:

for j=1:m+1
  dp(j)=0;
  for k=1:m+1
    if k~=j
      f=1;
      for l=1:m+1
        if l~=k & l~=j
          f=f*(c-t(l))/(t(j)-t(l));  
        end;
      end;
      dp(j)=dp(j)+f/(t(j)-t(k));
    end;
  end;
end;

return;
