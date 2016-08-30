function hopf=p_tohopf(funcs,point,freqs)
%% convert point to Hopf bifurcation point
% function hopf_point=p_tohopf(funcs,point {,freqs})
% INPUT:
%   funcs problem functions
%	point with stability information 
%   optional freqs: frequency to be excluded from consideration (or
%   used otherwise, depending on input point type)
% OUTPUT:
%	hopf_point uncorrected starting guess for hopf point
%
% If point is of type 'hoho' freqs can be string 'omega1' or 'omega2': then
% omega1 or omega2 will be selected as Hopf frequency, if freqs is value,
% it will be excluded
%
%
% $Id$
%
%%
if nargin<3
    freqs=[];
end

hopf.kind='hopf';
hopf.parameter=point.parameter;

switch point.kind
    case {'stst','fold','hopf'}
        [stability, eigv] = p_stabil_ndde(funcs, point);
        l1 = stability.l1;
        if isempty(l1)
            error('P_TOHOPF: point does not contain stability information l1!');
        end
        if strcmp(point.kind,'hopf')
            % remove known imaginary pair
            freqs=[abs(point.omega),freqs];
        end
        if ~isempty(freqs)
            for i=1:length(freqs)
                [i1,i2]=min(abs(real(l1))+abs(abs(imag(l1)-freqs(i)))); %#ok<ASGLU>
                l1(i2)=-1;
                [i1,i2]=min(abs(real(l1))+abs(abs(imag(l1)+freqs(i)))); %#ok<ASGLU>
                l1(i2)=-1;
            end
        end
        % look for non-real complex pair closest to imaginary axis
        idx = (imag(l1)>0);
        l1=l1(idx);
        eigv = eigv(:, idx);
        if isempty(l1)
            error('P_TOHOPF: no good pair of complex roots found.');
        end
        [i1,i2]=min(abs(real(l1))); %#ok<ASGLU>
        omega=imag(l1(i2));
        v = eigv(:,i2);
        x=point.x;
    case 'psol'
        x=sum(point.profile,2)/size(point.profile,2);
        omega=2*pi/point.period;
    otherwise
        try
            conversion=str2func([point.kind,'_tohopf']);
            hopf=conversion(funcs,point,freqs);
            return
        catch ME
            error('p_tohopf: point type %s not supported',point.kind);
        end
end

hopf.x=x;
if exist('v', 'var')
    % By now, v stores all the (interleaved) time profiles for the most
    % instable eigenvalue. But only the initial value of each profile
    % is need as the time evolution is fully determined by the eigenvalue.
    ndim = size(point.x, 1);
    hopf.v = v(1:ndim);
else
    D=ch_matrix(funcs,x,point.parameter,1i*omega);
    [E1,E2]=eig(D);
    [i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
    hopf.v=E1(:,i2);
end
hopf.omega=omega;
end

