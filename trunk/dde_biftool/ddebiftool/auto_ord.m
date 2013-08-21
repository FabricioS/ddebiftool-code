function [ITM1]=auto_ord(N,TM,N1,TM1)

% function [itm1]=auto_ord(n,tm,n1,tm1)
% INPUT:
%	n size of tm
%	tm ascending array 
%	n1 size of tm1
%	tm1 ascending array
% OUTPUT:
%	itm1 itm1(k) gives index of tm1(k) in tm-interval
% COMMENT:
%       this function is a matlab translation of the AUTO
%       fortran function ORDR 

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

K0=2;
for J1=1:N1,
  for J=K0:N,
    K1=J;
    if (TM1(J1)<TM(J)), break; end;
  end;
  ITM1(J1)=K1-1;
  K0=K1;
end;

return;
