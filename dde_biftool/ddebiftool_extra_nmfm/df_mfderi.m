function y = df_mfderi(funcs,xx,par,varargin)
% function y = df_mfderi(funcs,xx,par,varargin)
% PURPOSE: computes multilinear forms using finite differences
% INPUT:
%  funcs: structure containing (at least) sys_rhs and sys_deri
%  xx: state vectors
%  par: parameters
%  varargin: application vectors (1-5)
% OUTPUT:
%  y: result vector

if nargin == 3
	error('SYS_MFDERI: no arguments.');
elseif nargin > 8
	error('SYS_MFDERI: too many arguments.');
end

numarg = nargin - 3;

[n,r] = size(xx);
h = 10^(-6+numarg);

myrhs = funcs.sys_rhs;

switch numarg
	case 1
		u1 = varargin{1};
      y = 0*u1(:,1);
      for k = 0:r-1
         y = y + funcs.sys_deri(xx,par,k,[],[])*u1(:,k+1);
      end
	case 2
      h = 1e-4;
		u1 = varargin{1}; u2 = varargin{2};
		y = (B_eq(myrhs,h,u1+u2,xx,par) - B_eq(myrhs,h,u1-u2,xx,par))/4;
	case 3
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3};
      h = 1e-3/nthroot(norm(u1)*norm(u2)*norm(u3),3);
      %h = 1e-4;
		%y = (C_21(myrhs,h,u1+u3,u2+u3,xx,par) - C_21(myrhs,h,u1-u3,u2-u3,xx,par) - 2*C_eq(myrhs,h,u3,xx,par) - 2*C_21(myrhs,h,u1,u3,xx,par) )/4;
      y = (C_eq(myrhs,h,u1+u2+u3,xx,par) - C_eq(myrhs,h,u1+u2-u3,xx,par) ... 
         - C_eq(myrhs,h,u1-u2+u3,xx,par) + C_eq(myrhs,h,u1-u2-u3,xx,par))/24;
	case 4
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4};
      h = 1e-2/nthroot(norm(u1)*norm(u2)*norm(u3)*norm(u4),4)
      %h = 1e-2;
		y = (D_eq(myrhs,h,u1+u2+u3+u4,xx,par) - D_eq(myrhs,h,u1+u2+u3-u4,xx,par) ...
         + D_eq(myrhs,h,u1+u2-u3-u4,xx,par) - D_eq(myrhs,h,u1+u2-u3+u4,xx,par) ...
         - D_eq(myrhs,h,u1-u2+u3+u4,xx,par) + D_eq(myrhs,h,u1-u2+u3-u4,xx,par) ...
         - D_eq(myrhs,h,u1-u2-u3-u4,xx,par) + D_eq(myrhs,h,u1-u2-u3+u4,xx,par))/192;
	case 5
      h = 1e-3;
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4}; u5 = varargin{5};
      %h = 1e-2/nthroot(norm(u1)*norm(u2)*norm(u3)*norm(u4)*norm(u5),5)
		y = (E_eq(myrhs,h,u1+u2+u3+u4+u5,xx,par) - E_eq(myrhs,h,u1+u2+u3+u4-u5,xx,par) ...
         + E_eq(myrhs,h,u1+u2+u3-u4-u5,xx,par) - E_eq(myrhs,h,u1+u2+u3-u4+u5,xx,par) ...
         - E_eq(myrhs,h,u1+u2-u3+u4+u5,xx,par) + E_eq(myrhs,h,u1+u2-u3+u4-u5,xx,par) ...
         - E_eq(myrhs,h,u1+u2-u3-u4-u5,xx,par) + E_eq(myrhs,h,u1+u2-u3-u4+u5,xx,par) ...
         - E_eq(myrhs,h,u1-u2+u3+u4+u5,xx,par) + E_eq(myrhs,h,u1-u2+u3+u4-u5,xx,par) ...
         - E_eq(myrhs,h,u1-u2+u3-u4-u5,xx,par) + E_eq(myrhs,h,u1-u2+u3-u4+u5,xx,par) ...
         + E_eq(myrhs,h,u1-u2-u3+u4+u5,xx,par) - E_eq(myrhs,h,u1-u2-u3+u4-u5,xx,par) ...
         + E_eq(myrhs,h,u1-u2-u3-u4-u5,xx,par) - E_eq(myrhs,h,u1-u2-u3-u4+u5,xx,par))/1920;
end

   function z = B_eq(rhs,h,u,x0,p) % B(u,u)
      z = (rhs(x0 + h*u, p) + rhs(x0 - h*u, p))/h^2;
   end

   function z = C_eq(rhs,h,u,x0,p) % C(u,u,u)
      z = (rhs(x0 + 3*h*u, p) - 3*rhs(x0 + h*u, p) + 3*rhs(x0 - h*u,p) - rhs(x0 - 3*h*u,p))/(8*h^3);
   end

   function z = C_21(rhs,h,u,v,x0,p) % C(u,u,v)
      z = (C_eq(rhs,h,u+v,x0,p) - C_eq(rhs,h,u-v,x0,p))/6 - C_eq(rhs,h,v,x0,p)/3;
   end

   function z = D_eq(rhs,h,u,x0,p) % D(u,u,u,u)
      z = (rhs(x0 + 2*h*u,p) - 4*rhs(x0 + h*u,p) + 6*rhs(x0,p) - 4*rhs(x0 - h*u,p) + rhs(x0 - 2*h*u,p))/h^4;
   end

   function z = E_eq(rhs,h,u,x0,p) % E(u,u,u,u,u)
      z = (rhs(x0 + 5*h*u,p) - 5*rhs(x0 + 3*h*u,p) + 10*rhs(x0 + h*u,p) - 10*rhs(x0 - h*u,p) + 5*rhs(x0 - 3*h*u,p) - rhs(x0 - 5*h*u,p))/(32*h^5);
   end

end

