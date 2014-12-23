function D2Delta = nmfm_charmatderi2(funcs, xx,par,lambda)
%% Second derivative of characteristic matrix Delta(lambda) wrt lambda
% INPUT:
%   lambda: spectral value
% OUTPUT:
%   Delta: derivative of characteristic matrix at lambda
%
% $Id$
%
%%
[n,r] = size(xx); % n = #coordinates, r = #delays+1
eqtype = funcs.tp_del;
if eqtype == 0
   taus = par(sys_tau());
   taus = [0, taus]; % First delay zero
else % state-dependent delays
   taus = zeros(1,r); % First delay zero
   for i = 2:r
      taus(i) = sys_tau(i-1,xx,par);
   end
end
D2Delta = zeros(n);
for k = 1:r % For every delay
   Ak = funcs.sys_deri(xx,par,k-1,[],[]);
   D2Delta = D2Delta - taus(k)^2*Ak*exp(-lambda*taus(k));
end

end

