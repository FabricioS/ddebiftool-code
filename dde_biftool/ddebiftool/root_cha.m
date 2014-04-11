function [D,dD]=root_cha(funcs,x,par,l)
%% Characteristic matrix and its derivatives
% function [D,dD]=ch_matrx(funcs,x,par,l)
% INPUT:
%   funcs problem functions
%	x steady state solution in R^n
%	par parameter values
%	l root in C
% OUTPUT: 
%	D characteristic matrix in C^(n x n)
%	dD derivative of characteristic matrix wrt l in C^(n x n)

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% Id
%
%%
sys_tau=funcs.sys_tau;
sys_ntau=funcs.sys_ntau;
sys_deri=funcs.sys_deri;

tp_del=funcs.tp_del;
if tp_del==0
  tau=[0 par(sys_tau())];
  m=length(tau)-1;
  xx=x(:,ones(m+1,1));
else
  m=sys_ntau();
  xx=x(:,ones(m+1,1));
  t_tau=NaN(1,m);
  for j=1:m
    t_tau(j)=sys_tau(j,xx,par);
  end;
  tau=[0 t_tau];
end

n=length(x);

D=l*eye(n);
dD=eye(n);

for j=0:m
  B=sys_deri(xx,par,j,[],[])*exp(-l*tau(j+1));
  D=D-B;
  dD=dD+tau(j+1)*B;
end
end
