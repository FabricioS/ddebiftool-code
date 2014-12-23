function [p,q] = nmfm_border(A, p0, q0)
% PURPOSE: construct null vectors by bordering technique
% INPUT:
%	 A: n by n matrix
%   p0: approximate null vector s.t. p0*A = 0
%   q0: approximate null vector s.t. A*q0 = 0
% OUTPUT:
%   p: null vector s.t. p*A = 0
%   q: null vector s.t. A*q = 0

[n,~] = size(A);

if isempty(p0) || isempty(q0) % No previous null vectors specified
   Q = null(A); % Right
   P = null(A'); % Left
   if isempty(Q) || isempty(P)
      % Do it by eigenvalues, QZ method works best
      [V,D] = eig(A,eye(n),'qz');
      evs = diag(D);
      zeroeigind = find(abs(evs) < 1e-5); % the zero eigenvalue
      if length(zeroeigind) == 1
         Q = V(:,zeroeigind);
      end
      [V,D] = eig(A',eye(n),'qz');
      evs = diag(D);
      zeroeigind = find(abs(evs) < 1e-5); % the zero eigenvalue
      if length(zeroeigind) == 1
         P = V(:,zeroeigind);
      end
   end
   if isempty(Q) || isempty(P)
      % Matrix may be too close to singular, try rounding off to 4 decimals
      Ar = floor(1e4*A)*1e-4;
      Q = null(Ar);
      P = null(Ar');
      if isempty(Q) || isempty(P)
         fprintf('NMFM_BORDER: could not construct null vectors.\n');
         p = [];
         q = [];
         return;
      else
         fprintf('NMFM_BORDER: had to round off characteristic matrix.\n');
      end
   end
   pp0 = P(:,1)';
   qq0 = Q(:,1);
else
   pp0 = p0;
   qq0 = q0;
end
% Use bordering technique
Border = [A, pp0'; conj(qq0'), 0];
qsol = Border\[zeros(n,1) ;1];
q = qsol(1:n);
psol = [zeros(1,n), 1]/Border;
p = psol(1:n);


