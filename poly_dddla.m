function dp = poly_dddla(t,c)
% function dp=poly_dddla(t,c);
% INPUT:
%       t lagrange points in R^m+1
%       c evaluation point in R 
% OUTPUT:
%       dp values of third derivative of lagrange polynomials through t at c
%
% Copyright (c) David A.W. Barton (August 2005)

m = length(t) - 1;

for i = 1:(m+1)
    dp(i) = 0;
    for j = 1:(m + 1)
        if (j ~= i)
            dj = 0;
            for k = 1:(m + 1)
                if ((k ~= j) && (k ~= i))
                    df = 0;
                    for l = 1:(m + 1)
                        if ((l ~= k) && (l ~= j) && (l ~= i))
                            f = 1;
                            for r = 1:(m + 1)
                                if ((r ~= l) && (r ~= k) && (r ~= j) && (r ~= i))
                                    f = f*(c - t(r))/(t(i) - t(r));
                                end;
                            end;
                            df = df + f/(t(i) - t(l));
                        end;
                    end;
                    dj = dj + df/(t(i) - t(k));
                end;
            end;
            dp(i) = dp(i) + dj/(t(i) - t(j));
        end;
    end;
end;

return;
