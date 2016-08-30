function [stab,eigv] = p_stabil_ndde(funcs, st1, N, duration)
% function [stab,eigv] = p_stabil_ndde(funcs, st1, N)
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

sys_tau = funcs.sys_tau;
sys_deri = funcs.sys_deri;

if exist('N', 'var') ~= 1
    N = 100;
end;
if exist('duration', 'var') ~= 1
    duration = st1.parameter(sys_tau());
end
% dimension of system
n = size(st1.x,1);

% cheb returns the Chebyshev differentiation matrix [1..-1] so rescale to
% [-1..0]
M = -cheb(N)*2/duration;
A = zeros(n*(N+1));
for i=1:n
   A(i:n:end, i:n:end) = M; 
end

% Construct the solution vector
xx = st1.x*[1 1 0];

% Get the necessary derivatives
L0 = sys_deri(xx,st1.parameter,0,[],[]);
L1 = sys_deri(xx,st1.parameter,1,[],[]);
N1 = sys_deri(xx,st1.parameter,2,[],[]);

% Fill the last row block with the linearized problem
ilastrowblock = n*N+1:size(A, 1);
A(ilastrowblock, 1:n) = L1 + N1*M(1,1);
for j = 2:N
    A(ilastrowblock,(j-1)*n+(1:n)) = N1*M(1,j);
end;
A(ilastrowblock, n*N +(1:n)) = L0 + N1*M(1,end);

[eigv, l] = eig(A);
l = diag(l);
[~, idx] = sort(real(l), 1, 'descend');
eigv = eigv(:,idx);
stab.l0 = [];
stab.l1 = l(idx);
stab.n1 = [];
