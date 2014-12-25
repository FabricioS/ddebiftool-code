function hopf=hoho_tohopf(funcs,point,freqs) 
%% convert Hopf-Hopf point to Hopf bifurcation point
% function hopf_point=zeho_tohopf(funcs,point {,freqs})
% INPUT:
%   funcs problem functions
%	point with stability information 
%   optional freqs: frequency to be excluded/included
% OUTPUT:
%	hopf_point: Hopf point with single frequency
%
% $Id$
%
%%
hopf = struct(...
    'kind','hopf',...
    'parameter',point.parameter,...
    'x',point.x);
if isempty(freqs)
    omega = point.omega1;
elseif ischar(freqs)
    omega=point.(freqs);
else
    d1=min(abs(point.omega1-freqs));
    d2=min(abs(point.omega2-freqs));
    if d1<d2
        omega=point.omega2;
    else
        omega=point.omega1;
    end
end
hopf.omega=omega;
D=root_cha(funcs,point.x,point.parameter,1i*omega);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
hopf.v=E1(:,i2);
if isfield(point,'stability');
    hopf.stability=point.stability;
end
end
