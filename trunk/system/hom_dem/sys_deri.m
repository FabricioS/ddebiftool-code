function J=sys_deri(xx,par,nx,np,v)

% function J=sys_deri(xx,par,nx,np,v)
% INPUT:
%	xx state variable and delayed state variables columnwise
%	par list of parameter values
%	nx empty or list of requested state-derivatives (numbers of delay or zero) 
%	np empty or list of requested parameter-derivatives 
%	v matrix to multiply result with
% OUTPUT:
%	J result of derivatives on righthandside multiplied with v
% COMMENT:
%	the numerical derivatives are evaluated using forward differences

% (c) DDE-BIFTOOL v. 1.00, 11/03/2000

% first order derivative discretisation parameters:

abs_eps_x1=1e-6;
abs_eps_x2=1e-6;
abs_eps_p1=1e-6;
abs_eps_p2=1e-6;
rel_eps_x1=1e-6;
rel_eps_x2=1e-6;
rel_eps_p1=1e-6;
rel_eps_p2=1e-6;

n=size(xx,1);

J=[];

% first order derivatives of the state:
if length(nx)==1 & length(np)==0 & isempty(v),
  if nx(1)==0
    J(1,1)=-1;
    J(2,2)=-1;
  elseif nx(1)==1
    de=exp(-4*xx(1,2));
    dg=4*de/((1+de)^2);
    J(1,1)=par(1)*dg;
    J(1,2)=-par(2);
    J(2,1)=par(3)*dg;
  end;
% first order parameter derivatives:
elseif length(nx)==0 & length(np)==1 & isempty(v),
  f=sys_rhs(xx,par);
  par_eps=par;
  eps=abs_eps_p1+rel_eps_p1*abs(par(np));
  par_eps(np)=par(np)+eps;
  J=(sys_rhs(xx,par_eps)-f)/eps;
% second order state derivatives:
elseif length(nx)==2 & length(np)==0 & ~isempty(v),
  for j=1:n
    J(:,j)=sys_deri(xx,par,nx(1),[],[])*v;
    xx_eps=xx;
    eps=abs_eps_x2+rel_eps_x2*abs(xx(j,nx(2)+1));
    xx_eps(j,nx(2)+1)=xx_eps(j,nx(2)+1)+eps;
    J(:,j)=(sys_deri(xx_eps,par,nx(1),[],[])*v-J(:,j))/eps;
  end;
% mixed state parameter derivatives:
elseif length(nx)==1 & length(np)==1 & isempty(v),
  J=sys_deri(xx,par,nx(1),[],[]);
  par_eps=par;
  eps=abs_eps_p2+rel_eps_p2*abs(par(np));
  par_eps(np)=par(np)+eps;
  J=(sys_deri(xx,par_eps,nx(1),[],[])-J)/eps;
end;

if isempty(J)
  [nx np size(v)]
  error('SYS_DERI: requested derivative does not exist!');
end;

return;

