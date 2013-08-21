function f=sys_rhs(xx,par)

% q11 q12 q21 e1 e2 tau

g=1/(1+exp(-4*xx(1,2)));

f(1,1)=-xx(1,1)+par(1)*g-par(2)*xx(2,2)+par(4);
f(2,1)=-xx(2,1)+par(3)*g+par(5);

return;
