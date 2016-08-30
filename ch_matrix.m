function Delta=ch_matrix(funcs,xx,par,lambda,varargin)
%% Characteristic matrix and its derivatives
% function Delta=ch_matrx(funcs,x,par,l)
% INPUT:
%   funcs problem functions
%	x steady state solution in R^n (either n x (ntau+1), or n x 1)
%	par parameter values
%	lamba complex number at which  charactersitic matrix is computed
%   optional named argument: 'deri', integer (default 0) return derivative
%   of characteristic matrix wrt lambda
% OUTPUT: 
%	D characteristic matrix in C^(n x n)
%
% |xx| is assumed to be equilibrium such that |xx(:,2:end)| are ignored and
% |xx(:,ones(1,(ntau+1)) zeros(n, ntau)| is used instead.
%
% (c) DDE-BIFTOOL v. 2.00, 23/11/2001
%
% $Id$
%
%%
sys_tau = funcs.sys_tau;
sys_deri = funcs.sys_deri;
tp_del = funcs.tp_del;


default={'deri',0};
options=dde_set_options(default,varargin);
n = size(xx,1); % n = #coordinates,
if ~funcs.tp_del
    taus = par(sys_tau());
    taus = [0, taus]; % First delay zero
    r=length(taus); % number of delays, r = #delays+1
else % state-dependent delays
    r=funcs.sys_ntau()+1;
    taus = zeros(1,r); % First delay zero
    for i = 2:r
        taus(i) = sys_tau(i-1,xx(:,ones(1,i)),par);
    end
end

xx=[xx(:,ones(1,r)) zeros(n, r-1)];
nderi = options.deri;
lfac=[lambda, 1, zeros(1,nderi-1)];
Delta = lfac(nderi + 1)*eye(n);
if nderi==0
    Delta = Delta  - sys_deri(xx,par,0,[],[]);
end

for k = 2:r % For every delay (i.e. trivial loop by now)
    Ak = sys_deri(xx,par,k-1,[],[]);
    Bk = sys_deri(xx,par,k+r-2,[],[])*0;
    tauk = taus(k);
    % Successive derivatives of lambda*exp(-lambda*tauk) wrt lambda are
    % (-tauk)^(n-1)*(n-lambda*tauk)*exp(-lambda*tauk)
    Delta = Delta +(-tauk)^(nderi-1)*exp(-lambda*tauk)*(tauk*Ak + (nderi-lambda*tauk)*Bk);
end

end
