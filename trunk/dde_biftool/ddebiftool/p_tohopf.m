function hopf=p_tohopf(funcs,point,excludefreqs)
%% convert point to Hopf bifurcation point
% function hopf_point=p_tohopf(funcs,point {,excludefreqs})
% INPUT:
%   funcs problem functions
%	point with stability information 
% OUTPUT:
%	hopf_point uncorrected starting guess for hopf point

% (c) DDE-BIFTOOL v. 1.01, 14/07/2000
%
% $Id$
%
%%
if nargin<3
    excludefreqs=[];
end

hopf.kind='hopf';
hopf.parameter=point.parameter;

switch point.kind
    case {'stst','fold','hopf'}
        if ~isfield(point,'stability') || isempty(point.stability)
            error('P_TOHOPF: point does not contain stability information!');
        end
        l1=point.stability.l1;
        if isempty(l1)
            error('P_TOHOPF: point does not contain stability information l1!');
        end
        if strcmp(point.kind,'hopf')
            % remove known imaginary pair
            excludefreqs=[abs(point.omega),excludefreqs];
        end
        if ~isempty(excludefreqs)
            for i=1:length(excludefreqs)
                [i1,i2]=min(abs(real(l1))+abs(abs(imag(l1)-excludefreqs(i)))); %#ok<ASGLU>
                l1(i2)=-1;
                [i1,i2]=min(abs(real(l1))+abs(abs(imag(l1)+excludefreqs(i)))); %#ok<ASGLU>
                l1(i2)=-1;
            end
        end
        % look for non-real complex pair closest to imaginary axis
        l1=l1(imag(l1)>0);
        if isempty(l1)
            error('P_TOHOPF: no good pair of complex roots found.');
        end
        [i1,i2]=min(abs(real(l1))); %#ok<ASGLU>
        omega=imag(l1(i2));
        x=point.x;
    case 'psol'
        x=sum(point.profile,2)/size(point.profile,2);
        omega=2*pi/point.period;
    otherwise
        error('p_tohopf: point type %s not supported',point.kind);
end

hopf.x=x;
D=root_cha(funcs,x,point.parameter,1i*omega);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
hopf.v=E1(:,i2);
hopf.omega=omega;
end

