function [stab,eigv] = p_stabil_ndde(st1,N)
% function [stab,eigv] = p_stabil_ndde(st1,N)
%
% Use the method described in Breda, Maset, and Vermiglio 2006
% (Pseudospectral approximation of eigenvalues of derivative operators with
% non-local boundary conditions)
%
% st1 is a steady state solution
% N is the (optional) size of the discretisation to use
%
% Only works for single fixed delay (in current implementation)
% - general method works for multiple fixed delays

if exist('N') ~= 1
    N = 100;
end;

% dimension of system
n = size(st1.x,1);

% cheb returns the Chebyshev differentiation matrix [1..-1] so rescale to
% [-1..0]
M = -cheb(N)*2/st1.parameter(sys_tau());
% Fudge it to the correct dimension
N12 = (N+1)^2;
nN1 = n*(N+1);
M3 = [];
for i = 1:n
    M2 = [zeros(i-1,N12); reshape(M,1,N12); zeros(n-i,N12)];
    M3 = [M3; reshape(M2,nN1,N+1)];
end;
A = reshape(M3,nN1,nN1);

% Construct the solution vector
xx = st1.x*[1 1 0];

% Get the necessary derivatives
L0 = sys_deri(xx,st1.parameter,0,[],[]);
L1 = sys_deri(xx,st1.parameter,1,[],[]);
N1 = sys_deri(xx,st1.parameter,2,[],[]);

% Construct the A matrix
A(end-n+1:end,1:n) = L1 + N1*A(1,1);
for j = 2:N
    A(end-n+1:end,(j-1)*n+[1:n]) = N1*A(1,(j-1)*n+1);
end;
A(end-n+1:end,end-n+1:end) = L0 + N1*A(1,end-n+1);

if nargout == 2
    [eigv,l] = eig(A);
    l = diag(l);
    [ll,idx] = sort(real(l),1,'descend');
    eigv = eigv(:,idx);
else
    l = eig(A);
    [ll,idx] = sort(real(l),1,'descend');
end;

stab.l0 = [];
stab.l1 = l(idx);
stab.n1 = [];
