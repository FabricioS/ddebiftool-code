function [D,dD]=ch_matrx(x,par,l)

% function [D,dD]=ch_matrx(x,par,l)
% INPUT:
%	x steady state solution in R^n
%	par parameter values
%	l root in C
% OUTPUT: 
%	D characteristic matrix in C^(n x n)
%	dD derivative of characteristic matrix wrt l in C^(n x n)

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

tp_del=nargin('sys_tau');
if tp_del==0
  tau=[0 par(sys_tau)];
  m=length(tau)-1;
  xx=x;
  for j=1:m
    xx=[xx x];
  end;
else
  m=sys_ntau;
  xx=x;
  for j=1:m
    t_tau(j)=sys_tau(j,xx,par);
    xx=[xx x];
  end;
  tau=[0 t_tau];
end;

n=length(x);

D=l*eye(n);
dD=eye(n);

for j=0:m
  B=sys_deri(xx,par,j,[],[])*exp(-l*tau(j+1));
  D=D-B;
  dD=dD+tau(j+1)*B;
end;

return;


