function [J,res]=fold_jac(x,v,par,free_par,c)

% function [J,res]=fold_jac(x,v,par,free_par,c)
% INPUT:
%	x current fold solution guess in R^n
%	v current eigenvector guess in R^n
%	par current parameter values
%	free_par free parameter numbers
%	c normalization vector of v in R^(1 x n)
% OUTPUT: 
%	J jacobian in R^(2n+1+s x 2n+p)
%	res residual in R^(2n+1+s x 1)
 
% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

n=length(x);

tp_del=nargin('sys_tau');
if tp_del==0
  tau=par(sys_tau);
  m=length(tau);
  xx=x;
  for j=1:m
    xx=[xx x];
  end;
else
  m=sys_ntau;
  xx=x;
  for j=1:m
    tau(j)=sys_tau(j,xx,par);
    xx=[xx x];
  end;
end;

D=zeros(n);

for j=0:m
  B=sys_deri(xx,par,j,[],[]);
  D=D-B;
end;

Dv=D*v;
cv=c*v-1;
dDdxv=zeros(n);
for j=0:m
  for k=0:m % double work
    dDdxv=dDdxv-sys_deri(xx,par,[j k],[],v);
  end;
end;

dDdpv=zeros(n,length(free_par));

for k=1:length(free_par)
  for j=0:m
    dDdp=sys_deri(xx,par,j,free_par(k),[]);
    dDdpv(:,k)=dDdpv(:,k)-dDdp*v;
  end;
end;

res(1:n,1)=sys_rhs(xx,par);
res(n+1:2*n,1)=Dv;
res(2*n+1)=cv;

J=zeros(n,n);
for i=0:m
  J=J+sys_deri(xx,par,i,[],[]);
end;
for j=1:length(free_par)
  J(1:n,2*n+j)=sys_deri(xx,par,[],free_par(j),[]);
end;
J(n+1:2*n,1:n)=dDdxv;
J(n+1:2*n,n+1:2*n)=D;
J(n+1:2*n,2*n+(1:length(free_par)))=dDdpv;
J(2*n+1,n+1:2*n)=c;

return;

