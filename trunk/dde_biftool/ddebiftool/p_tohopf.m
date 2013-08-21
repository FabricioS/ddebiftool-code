function hopf=p_tohopf(point)

% function hopf_point=p_tohopf(point)
% INPUT:
%	point with stability information 
% OUTPUT:
%	hopf_point uncorrected starting guess for hopf point

% (c) DDE-BIFTOOL v. 1.01, 14/07/2000

hopf.kind='hopf';
hopf.parameter=point.parameter;

if point.kind=='stst' | point.kind=='fold' | point.kind=='hopf'
  if ~isfield(point,'stability')
    error('P_TOHOPF: point does not contain stability information!');
  end;
  if isempty(point.stability)
    error('P_TOHOPF: point does not contain stability information!');
  end;
  l1=point.stability.l1;
  if isempty(l1)
    error('P_TOHOPF: point does not contain enough stability information!');
  end;
  if point.kind=='hopf';
    % remove known imaginary pair
    [i1,i2]=min(abs(real(l1))+abs(abs(imag(l1))-abs(point.omega)));
    l1(i2)=0; 
    [i1,i2]=min(abs(real(l1))+abs(abs(imag(l1))-abs(point.omega)));
    l1(i2)=0;
  end;
  % look for non-real complex pair closest to imaginary axis
  m=max(abs(real(l1)));
  not_found=1;
  while not_found>0
    [i1,i2]=min(abs(real(l1)));
    omega=abs(imag(l1(i2)));
    if omega>0
      not_found=0;
    elseif abs(real(l1(i2)))>m+1
      not_found=-1;
    else
      l1(i2)=-m-2;
    end;
  end;
  x=point.x;
  if not_found==-1
    error('P_TOHOPF: no good pair of complex roots found.');
  end;
elseif point.kind=='psol'
  x=sum(point.profile,2)/size(point.profile,2);
  omega=2*pi/point.period;
end;

hopf.x=x;
[D,dD]=root_cha(x,point.parameter,sqrt(-1)*omega);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2)));
hopf.v=E1(:,i2);
hopf.omega=omega;

return;

