function Delta=ch_matrix(funcs,xx,par,lambda,varargin)
%% Characteristic matrix and its derivatives
% function Delta=ch_matrx(funcs,x,par,l)
% INPUT:
%   funcs problem functions
%	x steady state solution in R^n
%	par parameter values
%	l root in C
%   optional named argument: 'deri', integer (default 0) return derivative
%   of characteristic matrix wrt lambda
% OUTPUT: 
%	D characteristic matrix in C^(n x n)

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% $Id$
%
%%
default={'deri',0};
options=dde_set_options(default,varargin);
[n,r] = size(xx); % n = #coordinates, r = #delays+1
if ~funcs.tp_del
   taus = par(funcs.sys_tau());
   taus = [0, taus]; % First delay zero
else % state-dependent delays
   taus = zeros(1,r); % First delay zero
   for i = 2:r
      taus(i) = funcs.sys_tau(i-1,xx(:,1:i-1),par);
   end
end
lfac=[lambda,1,zeros(1,options.deri-1)];
tpow=options.deri;
tfac=(-1)^options.deri;
Delta = lfac(options.deri+1)*eye(n);
for k = 1:r % For every delay
   Ak = funcs.sys_deri(xx,par,k-1,[],[]);
   Delta = Delta -tfac*taus(k)^tpow*Ak*exp(-lambda*taus(k));
end

end
