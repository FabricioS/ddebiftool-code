function hopf=p_tohopf(funcs,point,freqs)
%% convert point to Hopf bifurcation point
% function hopf_point=p_tohopf(funcs,point {,freqs})
% INPUT:
%   funcs problem functions
%	point with stability information 
%   optional freqs: frequency to be excluded from consideration
% OUTPUT:
%	hopf_point uncorrected starting guess for hopf point
%
% If point is of type 'hoho' freqs can be string 'omega1' or 'omega2': then
% omega1 or omega2 will be selected as Hopf frequency, if freqs is value,
% it will be excluded

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
        if ~isfield(point,'stability') || isempty(point.stability)
            error('P_TOHOPF: point does not contain stability information!');
        end
        l1=point.stability.l1;
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
   % BW: Addition
   case {'genh','zeho'}
      hopf = point;
      x = point.x;
      omega = point.omega;
   % BW: Addition (extended by JS)
   case 'hoho'
      hopf = point;
      x = point.x;
      % this ensures hopf = p_tohopf(p_tohoho(hopf)):
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
      hopf = rmfield(hopf,'omega1');
      hopf = rmfield(hopf,'omega2');
      hopf.kind = 'hopf';        
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

