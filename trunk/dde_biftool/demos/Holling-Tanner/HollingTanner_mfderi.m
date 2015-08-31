function y = HollingTanner_mfderi(xx,par,varargin)
%% higher-order derivatives wrt to states (incl mixed)
%
% $Id$
%

if nargin == 2
	error('SYS_MFDERI: no arguments.');
elseif nargin > 7
	error('SYS_MFDERI: too many arguments.');
end

y = 0;

numarg = nargin - 2;

switch numarg
	case 1
		u1 = varargin{1};
		y = [(1-2*xx(1,1)-2*par(4)-xx(2,1)/(par(3)*xx(2,1)+xx(1,1))+xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2)*u1(1,1)+(-xx(1,1)/(par(3)*xx(2,1)+xx(1,1))+xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2*par(3))*u1(2,1);
			par(6)*(par(1)-xx(2,2)/xx(1,2))*u1(2,1)+par(6)*xx(2,1)*xx(2,2)/xx(1,2)^2*u1(1,2)-par(6)*xx(2,1)/xx(1,2)*u1(2,2)];
	case 2
		u1 = varargin{1}; u2 = varargin{2};
		y = [(-2+2*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2-2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3)*u1(1,1)*u2(1,1)+(-1/(par(3)*xx(2,1)+xx(1,1))+xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2*par(3)+xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^2-2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3))*u1(1,1)*u2(2,1)+(-1/(par(3)*xx(2,1)+xx(1,1))+xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^2*par(3)+xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^2-2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3))*u1(2,1)*u2(1,1)+(2*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^2*par(3)-2*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2)*u1(2,1)*u2(2,1);
			par(6)*xx(2,2)/xx(1,2)^2*u1(2,1)*u2(1,2)-par(6)/xx(1,2)*u1(2,1)*u2(2,2)+par(6)*xx(2,2)/xx(1,2)^2*u1(1,2)*u2(2,1)-par(6)/xx(1,2)*u1(2,2)*u2(2,1)-2*par(6)*xx(2,1)*xx(2,2)/xx(1,2)^3*u1(1,2)*u2(1,2)+par(6)*xx(2,1)/xx(1,2)^2*u1(1,2)*u2(2,2)+par(6)*xx(2,1)/xx(1,2)^2*u1(2,2)*u2(1,2)];
	case 3
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3};
		y = [(-6*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3+6*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4)*u1(1,1)*u2(1,1)*u3(1,1)+(2/(par(3)*xx(2,1)+xx(1,1))^2-4*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)-2*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^3+6*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3))*u1(1,1)*u2(1,1)*u3(2,1)+(2/(par(3)*xx(2,1)+xx(1,1))^2-4*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)-2*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^3+6*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3))*u1(1,1)*u2(2,1)*u3(1,1)+(2/(par(3)*xx(2,1)+xx(1,1))^2*par(3)-2*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2-4*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+6*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2)*u1(1,1)*u2(2,1)*u3(2,1)+(2/(par(3)*xx(2,1)+xx(1,1))^2-4*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)-2*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^3+6*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3))*u1(2,1)*u2(1,1)*u3(1,1)+(2/(par(3)*xx(2,1)+xx(1,1))^2*par(3)-2*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2-4*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+6*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2)*u1(2,1)*u2(1,1)*u3(2,1)+(2/(par(3)*xx(2,1)+xx(1,1))^2*par(3)-2*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2-4*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+6*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2)*u1(2,1)*u2(2,1)*u3(1,1)+(-6*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2+6*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3)*u1(2,1)*u2(2,1)*u3(2,1);
			-2*par(6)*xx(2,2)/xx(1,2)^3*u1(2,1)*u2(1,2)*u3(1,2)+par(6)/xx(1,2)^2*u1(2,1)*u2(1,2)*u3(2,2)+par(6)/xx(1,2)^2*u1(2,1)*u2(2,2)*u3(1,2)-2*par(6)*xx(2,2)/xx(1,2)^3*u1(1,2)*u2(2,1)*u3(1,2)+par(6)/xx(1,2)^2*u1(1,2)*u2(2,1)*u3(2,2)+par(6)/xx(1,2)^2*u1(2,2)*u2(2,1)*u3(1,2)-2*par(6)*xx(2,2)/xx(1,2)^3*u1(1,2)*u2(1,2)*u3(2,1)+par(6)/xx(1,2)^2*u1(1,2)*u2(2,2)*u3(2,1)+par(6)/xx(1,2)^2*u1(2,2)*u2(1,2)*u3(2,1)+6*par(6)*xx(2,1)*xx(2,2)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(1,2)-2*par(6)*xx(2,1)/xx(1,2)^3*u1(1,2)*u2(1,2)*u3(2,2)-2*par(6)*xx(2,1)/xx(1,2)^3*u1(1,2)*u2(2,2)*u3(1,2)-2*par(6)*xx(2,1)/xx(1,2)^3*u1(2,2)*u2(1,2)*u3(1,2)];
	case 4
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4};
		y = [(24*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5)*u1(1,1)*u2(1,1)*u3(1,1)*u4(1,1)+(-6/(par(3)*xx(2,1)+xx(1,1))^3+18*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)+6*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3))*u1(1,1)*u2(1,1)*u3(1,1)*u4(2,1)+(-6/(par(3)*xx(2,1)+xx(1,1))^3+18*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)+6*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3))*u1(1,1)*u2(1,1)*u3(2,1)*u4(1,1)+(-8/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+12*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2+12*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2)*u1(1,1)*u2(1,1)*u3(2,1)*u4(2,1)+(-6/(par(3)*xx(2,1)+xx(1,1))^3+18*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)+6*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3))*u1(1,1)*u2(2,1)*u3(1,1)*u4(1,1)+(-8/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+12*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2+12*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2)*u1(1,1)*u2(2,1)*u3(1,1)*u4(2,1)+(-8/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+12*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2+12*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2)*u1(1,1)*u2(2,1)*u3(2,1)*u4(1,1)+(-6/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2+6*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3+18*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3)*u1(1,1)*u2(2,1)*u3(2,1)*u4(2,1)+(-6/(par(3)*xx(2,1)+xx(1,1))^3+18*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)+6*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3))*u1(2,1)*u2(1,1)*u3(1,1)*u4(1,1)+(-8/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+12*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2+12*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2)*u1(2,1)*u2(1,1)*u3(1,1)*u4(2,1)+(-8/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+12*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2+12*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2)*u1(2,1)*u2(1,1)*u3(2,1)*u4(1,1)+(-6/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2+6*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3+18*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3)*u1(2,1)*u2(1,1)*u3(2,1)*u4(2,1)+(-8/(par(3)*xx(2,1)+xx(1,1))^3*par(3)+12*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2+12*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2)*u1(2,1)*u2(2,1)*u3(1,1)*u4(1,1)+(-6/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2+6*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3+18*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3)*u1(2,1)*u2(2,1)*u3(1,1)*u4(2,1)+(-6/(par(3)*xx(2,1)+xx(1,1))^3*par(3)^2+6*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3+18*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3)*u1(2,1)*u2(2,1)*u3(2,1)*u4(1,1)+(24*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3-24*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^4)*u1(2,1)*u2(2,1)*u3(2,1)*u4(2,1);
			6*par(6)*xx(2,2)/xx(1,2)^4*u1(2,1)*u2(1,2)*u3(1,2)*u4(1,2)-2*par(6)/xx(1,2)^3*u1(2,1)*u2(1,2)*u3(1,2)*u4(2,2)-2*par(6)/xx(1,2)^3*u1(2,1)*u2(1,2)*u3(2,2)*u4(1,2)-2*par(6)/xx(1,2)^3*u1(2,1)*u2(2,2)*u3(1,2)*u4(1,2)+6*par(6)*xx(2,2)/xx(1,2)^4*u1(1,2)*u2(2,1)*u3(1,2)*u4(1,2)-2*par(6)/xx(1,2)^3*u1(1,2)*u2(2,1)*u3(1,2)*u4(2,2)-2*par(6)/xx(1,2)^3*u1(1,2)*u2(2,1)*u3(2,2)*u4(1,2)-2*par(6)/xx(1,2)^3*u1(2,2)*u2(2,1)*u3(1,2)*u4(1,2)+6*par(6)*xx(2,2)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(2,1)*u4(1,2)-2*par(6)/xx(1,2)^3*u1(1,2)*u2(1,2)*u3(2,1)*u4(2,2)-2*par(6)/xx(1,2)^3*u1(1,2)*u2(2,2)*u3(2,1)*u4(1,2)-2*par(6)/xx(1,2)^3*u1(2,2)*u2(1,2)*u3(2,1)*u4(1,2)+6*par(6)*xx(2,2)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(1,2)*u4(2,1)-2*par(6)/xx(1,2)^3*u1(1,2)*u2(1,2)*u3(2,2)*u4(2,1)-2*par(6)/xx(1,2)^3*u1(1,2)*u2(2,2)*u3(1,2)*u4(2,1)-2*par(6)/xx(1,2)^3*u1(2,2)*u2(1,2)*u3(1,2)*u4(2,1)-24*par(6)*xx(2,1)*xx(2,2)/xx(1,2)^5*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)+6*par(6)*xx(2,1)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(1,2)*u4(2,2)+6*par(6)*xx(2,1)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(2,2)*u4(1,2)+6*par(6)*xx(2,1)/xx(1,2)^4*u1(1,2)*u2(2,2)*u3(1,2)*u4(1,2)+6*par(6)*xx(2,1)/xx(1,2)^4*u1(2,2)*u2(1,2)*u3(1,2)*u4(1,2)];
	case 5
		u1 = varargin{1}; u2 = varargin{2}; u3 = varargin{3}; u4 = varargin{4}; u5 = varargin{5};
		y = [(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(2,1)*u2(1,1)*u3(2,1)*u4(1,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(2,1)*u2(1,1)*u3(2,1)*u4(1,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(2,1)*u2(1,1)*u3(2,1)*u4(2,1)*u5(1,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3-24*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^4-96*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^4)*u1(2,1)*u2(1,1)*u3(2,1)*u4(2,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(2,1)*u2(2,1)*u3(1,1)*u4(1,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(2,1)*u2(2,1)*u3(1,1)*u4(1,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(2,1)*u2(2,1)*u3(1,1)*u4(2,1)*u5(1,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3-24*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^4-96*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^4)*u1(2,1)*u2(2,1)*u3(1,1)*u4(2,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(2,1)*u2(2,1)*u3(2,1)*u4(1,1)*u5(1,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3-24*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^4-96*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^4)*u1(2,1)*u2(2,1)*u3(2,1)*u4(1,1)*u5(2,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3-24*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^4-96*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^4)*u1(2,1)*u2(2,1)*u3(2,1)*u4(2,1)*u5(1,1)+(-120*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^4+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^5)*u1(2,1)*u2(2,1)*u3(2,1)*u4(2,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(1,1)*u2(1,1)*u3(1,1)*u4(2,1)*u5(2,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4-96*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)-24*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3))*u1(1,1)*u2(1,1)*u3(2,1)*u4(1,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(1,1)*u2(1,1)*u3(2,1)*u4(1,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(1,1)*u2(1,1)*u3(2,1)*u4(2,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(1,1)*u2(1,1)*u3(2,1)*u4(2,1)*u5(2,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4-96*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)-24*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3))*u1(1,1)*u2(2,1)*u3(1,1)*u4(1,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(1,1)*u2(2,1)*u3(1,1)*u4(1,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(1,1)*u2(2,1)*u3(1,1)*u4(2,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(1,1)*u2(2,1)*u3(1,1)*u4(2,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(1,1)*u2(2,1)*u3(2,1)*u4(1,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(1,1)*u2(2,1)*u3(2,1)*u4(1,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(1,1)*u2(2,1)*u3(2,1)*u4(2,1)*u5(1,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^3-24*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^4-96*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^4)*u1(1,1)*u2(2,1)*u3(2,1)*u4(2,1)*u5(2,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4-96*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)-24*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3))*u1(2,1)*u2(1,1)*u3(1,1)*u4(1,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(2,1)*u2(1,1)*u3(1,1)*u4(1,1)*u5(2,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)-72*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2-48*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^2)*u1(2,1)*u2(1,1)*u3(1,1)*u4(2,1)*u5(1,1)+(36/(par(3)*xx(2,1)+xx(1,1))^4*par(3)^2-48*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^3-72*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)^2+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3)^3)*u1(2,1)*u2(1,1)*u3(1,1)*u4(2,1)*u5(2,1)+(-120*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6)*u1(1,1)*u2(1,1)*u3(1,1)*u4(1,1)*u5(1,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4-96*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)-24*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3))*u1(1,1)*u2(1,1)*u3(1,1)*u4(1,1)*u5(2,1)+(24/(par(3)*xx(2,1)+xx(1,1))^4-96*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^5*par(3)-24*xx(1,1)/(par(3)*xx(2,1)+xx(1,1))^5+120*xx(1,1)*xx(2,1)/(par(3)*xx(2,1)+xx(1,1))^6*par(3))*u1(1,1)*u2(1,1)*u3(1,1)*u4(2,1)*u5(1,1);
			120*par(6)*xx(2,1)*xx(2,2)/xx(1,2)^6*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(2,2)*u3(1,2)*u4(1,2)*u5(2,1)+6*par(6)/xx(1,2)^4*u1(2,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(2,1)+6*par(6)/xx(1,2)^4*u1(2,1)*u2(2,2)*u3(1,2)*u4(1,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(2,1)*u3(1,2)*u4(1,2)*u5(2,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(2,1)*u3(1,2)*u4(2,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(2,1)*u3(2,2)*u4(1,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(2,2)*u2(2,1)*u3(1,2)*u4(1,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(2,1)*u4(1,2)*u5(2,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(2,1)*u4(2,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(2,2)*u3(2,1)*u4(1,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(2,2)*u2(1,2)*u3(2,1)*u4(1,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(1,2)*u4(2,1)*u5(2,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(2,2)*u4(2,1)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(2,2)*u3(1,2)*u4(2,1)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(2,2)*u2(1,2)*u3(1,2)*u4(2,1)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(1,2)*u4(2,2)*u5(2,1)+6*par(6)/xx(1,2)^4*u1(1,2)*u2(1,2)*u3(2,2)*u4(1,2)*u5(2,1)+6*par(6)/xx(1,2)^4*u1(2,1)*u2(1,2)*u3(1,2)*u4(1,2)*u5(2,2)+6*par(6)/xx(1,2)^4*u1(2,1)*u2(1,2)*u3(1,2)*u4(2,2)*u5(1,2)+6*par(6)/xx(1,2)^4*u1(2,1)*u2(1,2)*u3(2,2)*u4(1,2)*u5(1,2)-24*par(6)*xx(2,2)/xx(1,2)^5*u1(2,1)*u2(1,2)*u3(1,2)*u4(1,2)*u5(1,2)-24*par(6)*xx(2,2)/xx(1,2)^5*u1(1,2)*u2(2,1)*u3(1,2)*u4(1,2)*u5(1,2)-24*par(6)*xx(2,2)/xx(1,2)^5*u1(1,2)*u2(1,2)*u3(2,1)*u4(1,2)*u5(1,2)-24*par(6)*xx(2,2)/xx(1,2)^5*u1(1,2)*u2(1,2)*u3(1,2)*u4(2,1)*u5(1,2)-24*par(6)*xx(2,2)/xx(1,2)^5*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(2,1)-24*par(6)*xx(2,1)/xx(1,2)^5*u1(1,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(2,2)-24*par(6)*xx(2,1)/xx(1,2)^5*u1(1,2)*u2(1,2)*u3(1,2)*u4(2,2)*u5(1,2)-24*par(6)*xx(2,1)/xx(1,2)^5*u1(1,2)*u2(1,2)*u3(2,2)*u4(1,2)*u5(1,2)-24*par(6)*xx(2,1)/xx(1,2)^5*u1(1,2)*u2(2,2)*u3(1,2)*u4(1,2)*u5(1,2)-24*par(6)*xx(2,1)/xx(1,2)^5*u1(2,2)*u2(1,2)*u3(1,2)*u4(1,2)*u5(1,2)];
end