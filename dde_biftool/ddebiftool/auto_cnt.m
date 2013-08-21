function [D]=auto_cnt(N)

% function [d]=auto_cnt(n)
% INPUT:
%	n n-th derivative
% OUTPUT:
%	d coefficients of central difference formula
% COMMENT: 
%       this function is a matlab translation of the AUTO
%       fortran function CNTDIF

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

D(1)=1;
       
if (N==0),
  return;
end;

for I=1:N,
  D(I+1)=0;
  for K=1:I,
    K1=I+2-K;
    D(K1)=D(K1-1)-D(K1);
  end;
  D(1)=-D(1);
end;

% Scale to [0,1]  :

SC=N^N;
NP1=N+1;
for I=1:NP1,
  D(I)=SC*D(I);
end;

return;
