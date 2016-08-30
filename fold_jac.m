function [J,res]=fold_jac(funcs,x,v,par,free_par,c)
%% Jacobian and residual of nonlinear system for fold
% function [J,res]=fold_jac(funcs,x,v,par,free_par,c)
% INPUT:
%   funcs problem functions
%   x current fold solution guess in R^n
%   v current eigenvector guess in R^n
%   par current parameter values
%   free_par free parameter numbers
%   c normalization vector of v in R^(1 x n)
% OUTPUT: 
%   J jacobian in R^(2n+1+s x 2n+p)
%   res residual in R^(2n+1+s x 1)
 
% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

%%
n=length(x);
sys_tau=funcs.sys_tau;
sys_rhs=funcs.sys_rhs;
sys_ntau=funcs.sys_ntau;
sys_deri=funcs.sys_deri;

tp_del=funcs.tp_del;
if tp_del==0
  tau = par(sys_tau());
  m=length(tau);
  xx = [x(:,ones(m+1,1)) zeros(n,m)];
else
    m=sys_ntau();
    xx=x(:,ones(m+1,1));
    tau=NaN(1,m);
    for j=1:m
        tau(j)=sys_tau(j,xx,par);
    end
end

n_fp=length(free_par);

% D: characteristic matrix
D=zeros(n);

% Derivative of D wrt to the equilibrium state
% dDxv.dx ~ (D([x, x...x, 0...0], par, l)-D([x-dx, x-dx...x-dx, 0...0], par, l)).v
dDdxv=zeros(n);

% Derivative of D wrt to free parameters
% dDdpv.[dpar_1...dpar_nfp] = (D([x,0], par, l) - D([x,0], par-dpar, l)).v
dDdpv=zeros(n,n_fp);

% Jacobian for the equilibrium state
Jx = zeros(n, n);
Jp = zeros(n, n_fp);

L0 = sys_deri(xx, par, 0, [], []);
Jx(:,:) = L0;
D(:,:) = - L0;
for k=0:m
  dDdxv = dDdxv - sys_deri(xx, par, [0 k], [], v);
end
for k=1:n_fp
  Jp(:,k) = sys_deri(xx,par,[],free_par(k),[]);
  dDdpv(:,k) = dDdpv(:,k)- sys_deri(xx,par,0,free_par(k),[]) * v;
end

for j=1:m
  Lj = sys_deri(xx,par,j,[],[]);
  Jx = Jx + Lj;
  D = D - Lj;

  for k=0:m
    dDdxv = dDdxv - sys_deri(xx,par,[j k],[],v);
  end;

  for k=1:n_fp
    dDdpv(:,k) = dDdpv(:,k) - sys_deri(xx,par,j,free_par(k),[]) * v;
  end;
end;

Dv=D*v;
cv=c*v-1;

res = [sys_rhs(xx,par);
       Dv;
       cv];

tmp1 = zeros(n,1);
tmpn = zeros(n,n);
tmpnf = zeros(n_fp, 1);
J = [Jx    , tmpn, Jp;
     dDdxv , D   , dDdpv;
     tmp1.', c   , tmpnf.'];

end
