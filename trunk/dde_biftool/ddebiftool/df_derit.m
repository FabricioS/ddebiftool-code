function dtau=sys_dtau(delay_nr,xx,par,nx,np)

% function dtau=sys_dtau(delay_nr,xx,par,nx,np)
% INPUT:
%       delay_nr delay number 
%	xx state variable and delayed state variables columnwise
%	par list of parameter values
%	nx empty or list of requested state-derivatives (numbers of delay or zero) 
%	np empty or list of requested parameter-derivatives 
% OUTPUT:
%	dtau result of derivatives on delay function
% COMMENT:
%	the numerical derivatives are evaluated using forward differences

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

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

dtau=[];

% first order state derivatives:
if length(nx)==1 & length(np)==0,
  tau=sys_tau(delay_nr,xx,par);
  for j=1:n
    xx_eps=xx;
    eps=abs_eps_x1+rel_eps_x1*abs(xx(j,nx+1));
    xx_eps(j,nx+1)=xx(j,nx+1)+eps;
    dtau(j)=(sys_tau(delay_nr,xx_eps,par)-tau)/eps;
  end;
% first order parameter derivatives:
elseif length(nx)==0 & length(np)==1,
  tau=sys_tau(delay_nr,xx,par);
  par_eps=par;
  eps=abs_eps_p1+rel_eps_p1*abs(par(np));
  par_eps(np)=par(np)+eps;
  dtau=(sys_tau(delay_nr,xx,par_eps)-tau)/eps;
% second order state derivatives:
elseif length(nx)==2 & length(np)==0,
    dtau_init=(sys_dtau(delay_nr,xx,par,nx(1),[]))';
  for j=1:n
    xx_eps=xx;
    eps=abs_eps_x2+rel_eps_x2*abs(xx(j,nx(2)+1));
    xx_eps(j,nx(2)+1)=xx_eps(j,nx(2)+1)+eps;
    dtau(:,j)=((sys_dtau(delay_nr,xx_eps,par,nx(1),[]))'-dtau_init)/eps;
  end;
% mixed state parameter derivatives:
elseif length(nx)==1 & length(np)==1,
  dtau=sys_dtau(delay_nr,xx,par,nx(1),[]);
  par_eps=par;
  eps=abs_eps_p2+rel_eps_p2*abs(par(np));
  par_eps(np)=par(np)+eps;
  dtau=(sys_dtau(delay_nr,xx,par_eps,nx(1),[])-dtau)/eps;
end;

if isempty(dtau)
  [delay_nr nx np]
  error('SYS_DTAU: requested derivative does not exist!');
end;

return;

