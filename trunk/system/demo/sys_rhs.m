function f=sys_rhs(xx,par)

% kappa beta a12 a21 tau1 tau2 tau_s

f(1,1)=-par(1)*xx(1,1)+par(2)*tanh(xx(1,4))+par(3)*tanh(xx(2,3));
f(2,1)=-par(1)*xx(2,1)+par(2)*tanh(xx(2,4))+par(4)*tanh(xx(1,2));

return;
