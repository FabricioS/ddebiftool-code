function dp=poly_ddel(m,c);

% function dp=poly_ddel(m,c);
% INPUT:
%       m lagrange degree
%       c evaluation point in [0,1]
% OUTPUT:
%       dp values of third derivative of lagrange polynomials through 0:1/m:1 at c

% (c) David Barton 19/11/2005
% Adapted from
% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

if m==1
  dp(1)=0;
  dp(2)=0;
elseif m==2
  dp(1)=0;
  dp(2)=0;
  dp(3)=0;
elseif m==3
  dp(1)=-27;
  dp(2)=81;
  dp(3)=-81;
  dp(4)=27;
elseif m==4
  dp(1)=256*c-160;
  dp(2)=-1024*c+576;
  dp(3)=1536*c-768;
  dp(4)=-1024*c+448;
  dp(5)=256*c-96;
else
  t = [0:m]/m;
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
end;

return;
