function PHI = nmfm_handletomatrix(fn, arg)
% INPUT:
%   fn: function handle
%   arg: argument vector for function
% OUTPUT:
%   PHI: n by r matrix

n = length(fn(0));
r = length(arg);

PHI = zeros(n,r);

for j = 1:n
    for k = 1:r
        phivec = fn(arg(k));
        PHI(j,k) = phivec(j);
    end
end


end

