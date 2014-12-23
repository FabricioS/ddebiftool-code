function DDelta = nmfm_charmatderi(funcs, xx,par,lambda)
% INPUT:
%   lambda: spectral value
% OUTPUT:
%   Delta: derivative of characteristic matrix at lambda

[n,r] = size(xx); % n = #coordinates, r = #delays+1
eqtype = nargin(funcs.sys_tau);
if eqtype == 0
   taus = par(funcs.sys_tau());
   taus = [0, taus]; % First delay zero
else % state-dependent delays
   taus = zeros(1,r); % First delay zero
   for i = 2:r
      taus(i) = sys_tau(i-1,xx,par);
   end
end

DDelta = eye(n);
for k = 1:r % For every delay
   Ak = funcs.sys_deri(xx,par,k-1,[],[]);
   DDelta = DDelta + taus(k)*Ak*exp(-lambda*taus(k));
end

end

