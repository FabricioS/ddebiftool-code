function [J,res]=stst_jac(x,par,free_par)

% function [J,res]=stst_jac(x,par,free_par)
% INPUT:
%	x current solution guess in R^n
%	par current parameter values
%	free_par free parameter numbers
% OUTPUT: 
%	J jacobian in R^(n+s x n+p)
%	res residual in R^(n+s x 1)

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

n=length(x); % system dimension 

tp_del=nargin('sys_tau');
if tp_del==0,
  m=length(sys_tau); % number of delays
else
  m=sys_ntau;
end;

xx=x;
for i=1:m
  xx=[xx x];
end;

res=sys_rhs(xx,par);

J=zeros(n,n+length(free_par));

for i=0:m
  J(1:n,1:n)=J(1:n,1:n)+sys_deri(xx,par,i,[],[]);
end;
for j=1:length(free_par)
  J(1:n,n+j)=J(1:n,n+j)+sys_deri(xx,par,[],free_par(j),[]);
end;

return;

