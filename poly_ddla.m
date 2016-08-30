function dp = poly_ddla(t,c)
% function dp=poly_ddla(t,c);
% INPUT:
%       t lagrange points in R^m+1
%       c evaluation point in R 
% OUTPUT:
%       dp values of second derivative of lagrange polynomials through t at c
%
% Copyright (c) David A.W. Barton (August 2005)

m = length(t) - 1;

for j = 1:(m + 1)
    dp(j) = 0;
    for k = 1:(m + 1)
        if (k ~= j)
            df = 0;
            for l = 1:(m + 1)
                if ((l ~= k) && (l ~= j))
                    f = 1;
                    for r = 1:(m + 1)
                        if ((r ~= l) && (r ~= k) && (r ~= j))
                            f = f*(c - t(r))/(t(j) - t(r));
                        end;
                    end;
                    df = df + f/(t(j) - t(l));
                end;
            end;
            dp(j) = dp(j) + df/(t(j) - t(k));
        end;
    end;
end;

return;
