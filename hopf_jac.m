function [J,res]=hopf_jac(funcs,x,omega,v,par,free_par,c)
%% Jacobian and residual of nonlinear system for Hopf
% function [J,res]=hopf_jac(funcs,x,omega,v,par,free_par,c)
% INPUT:
%   funcs problem function
%   x current Hopf solution guess in R^n
%   omega current Hopf frequency guess in R
%   v current eigenvector guess in C^n
%   par current parameters values in R^p
%   free_par free parameter numbers in N^d 
%   c normalization vector of v in C^(1 x n)
%OUTPUT: 
%   J jacobian in R^(3n+2+s x 3n+1+p)
%   res residual in R^(3n+2+s x 1)

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

% FIXME: handle delays that depends on neutral terms

n=length(x);
sys_tau=funcs.sys_tau;
sys_rhs=funcs.sys_rhs;
sys_ntau=funcs.sys_ntau;
sys_deri=funcs.sys_deri;
sys_dtau=funcs.sys_dtau;

tp_del=funcs.tp_del;
if tp_del==0
  n_tau = sys_tau();
  tau = par(n_tau);
  m=length(n_tau);
  xx = [x(:,ones(m+1,1)) zeros(n,m)];
else
  m=sys_ntau();
  xx = [x(:,ones(m+1,1)) zeros(n,m)];
  tau = NaN(1,m);
  for j=1:m
    tau(j)=sys_tau(j,xx,par);
  end
end

n_fp=length(free_par);
l=1i*omega;

% D: characteristic matrix
% dD: its derivative wrt to lambda(=1i*omega)
D = zeros(n);
dD=eye(n);

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
D(:,:) = l*dD - L0;
for k=0:m
  dDdxv = dDdxv-sys_deri(xx,par,[0 k],[],v);
end
for k=1:n_fp
  Jp(:,k) = sys_deri(xx,par,[],free_par(k),[]);
  dDdpv(:,k) = dDdpv(:,k)- sys_deri(xx,par,0,free_par(k),[]) * v;
end


for j=1:m
  expterm = exp(-l*tau(j));

  Lj = sys_deri(xx,par,j,[],[]);
  % Linearized neutral terms
  Nj = sys_deri(xx,par,m+j,[],[]);

  Jx = Jx + Lj;
  D=D - (Lj + l*Nj) * expterm;
  dD=dD + (tau(j)*Lj - (1-l*tau(j))*Nj) * expterm;

  for k=0:m
    dLjdNjxk = sys_deri(xx,par,[j k],[],v) + l*sys_deri(xx,par,[m+j k],[],v);
    dDdxv = dDdxv - dLjdNjxk * expterm;
    if tp_del~=0
      vdtau = v*sys_dtau(j,xx,par,k,[]);
      dDdxv = dDdxv + (Lj + l*Nj) * l * expterm * vdtau;
    end;
  end;

  for k=1:n_fp
    dLjdNjp = sys_deri(xx,par,j,free_par(k),[]) + l*sys_deri(xx,par,m+j,free_par(k),[]);
    dDdpv(:,k) = dDdpv(:,k) - dLjdNjp * v * expterm;
    if tp_del==0
      if n_tau(j)==free_par(k)
        dDdpv(:,k) = dDdpv(:,k) + (Lj + l*Nj) * v * l * expterm;
      end;
    elseif tp_del~=0
      vdtau = v*sys_dtau(j,xx,par,[],free_par(k)); % d=d(tau(j))/d(free_par(k))
      dDdpv(:,k)=dDdpv(:,k) + (Lj + l*Nj) * l * expterm * vdtau;
    end;
  end;
end;

Dv=D*v;
dDv = dD*v;
cv=c*v-1;

res = [sys_rhs(xx,par);
       real(Dv);
       imag(Dv);
       real(cv);
       imag(cv)];

tmp1 = zeros(n,1);
tmpn = zeros(n,n);
tmpnf = zeros(n_fp, 1);
J = [Jx         , tmpn   ,  tmpn   ,  tmp1     , Jp;
     real(dDdxv), real(D), -imag(D), -imag(dDv), real(dDdpv);
     imag(dDdxv), imag(D),  real(D),  real(dDv), imag(dDdpv);
     tmp1.'     , real(c), -imag(c),  0        , tmpnf.';
     tmp1.'     , imag(c),  real(c),  0        , tmpnf.'];

return;
