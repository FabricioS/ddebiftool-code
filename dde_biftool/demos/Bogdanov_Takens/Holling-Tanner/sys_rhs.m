function f = sys_rhs(xx,par)

f(1,1) = (xx(1,1)+par(4))*(1-xx(1,1)-par(4))-xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))-par(5);
f(2,1) = par(6)*xx(2,1)*(par(1)-xx(2,2)/xx(1,2));

end
