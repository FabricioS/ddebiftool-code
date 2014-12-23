function PHI = nmfm_handletomatrix(fn, arg)
%% Call single-argument function with array of arguments
% INPUT:
%   fn: function handle
%   arg: argument vector for function
% OUTPUT:
%   PHI: n by r matrix
%
% $Id$
%
%%
for k = length(arg):-1:1
    phivec = fn(arg(k));
    PHI(:,k) = phivec(:);
end
end

