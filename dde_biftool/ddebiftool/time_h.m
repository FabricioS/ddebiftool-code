function h=time_h(x,alpha,beta,real_part,h_min,h_max,par,rho_safety)

% function h=time_h(x,alpha,beta,real_part,h_min,h_max,par,rho_safety)
% INPUT:
%       x steady state solution in R^n 
%	alpha alpha-LMS parameters in R^k
%	beta beta-LMS parameters in R^k
%	real_part to compute roots with real part >= real_part  
%	h_min minimal h (relative to max(tau))
%	h_max maximal h (relative to max(tau))
%	par current parameter values in R^p
%	rho_safety safety radius of LMS-method
% OUTPUT: 
%       h recommended steplength (relative to max(tau))

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

tp_del=nargin('sys_tau');
if tp_del==0 % DDE case:
  tau=par(sys_tau);
  m=length(tau);
  xx=x;
  for i=1:m;
    xx=[xx x];
  end;
else % sd-DDE case:
  m=sys_ntau;
  xx=x;
  for j=1:m;
    tau(j)=sys_tau(j,xx,par);
    xx=[xx x];
  end;
end;

taumax=max(tau);
k=length(alpha);

D=sys_deri(xx,par,0,[],[]);
n0=norm(D);
normAB=n0+abs(real_part);
for i=1:m
  D=sys_deri(xx,par,i,[],[]);
  normAB=normAB+norm(D)*exp(-real_part*tau(i));
end;

ode=0;

% h1 for stability:

if normAB>0
  h1=0.9*rho_safety/normAB;
elseif normAB==0
  ode=1;
else % if normAB == inf
  h1=h_min*taumax/2;
end;

% h2 for invertibility:

if n0>0 
  h2=0.9*alpha(k)/(beta(k)*n0);
elseif n0==0
  h2=h1;
end;

% h:

if ~ode
  if taumax>0
    h=min(min([h1 h2])/taumax,h_max);
  else
    ode=1;
  end;
end;

if ode  
  disp('TIME_H warning: dde has become ode.');
  h=0;
elseif h<h_min
  h=h_min;
  if h_max>h_min
    disp('TIME_H warning: h_min is reached.');
  elseif h_max<h_min
    err=[h_min h_max]
    error('TIME_H: h_min>h_max not allowed.'); 
  end;
end;

return;

