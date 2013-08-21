function [l1,n1]=root_nwt(x,l0,max_n,epsi,par)

% function [l1,n]=root_nwt(x,l0,max_n,epsi,par)
% INPUT:
%	x steady state in RR^n
%	l0 number of root guesses in C^L
%	max_n maximum number of corrections
%	epsi stop criterion for newton iterations > 0
%       par current parameter values in R^p
% OUTPUT: 
%	l1 corrected roots in C^L
%	n1 number of corrections done (-1 if no convergence)

% (c) DDE-BIFTOOL v. 2.00, 23/11/2001

n=length(x);

tp_del=nargin('sys_tau');
if tp_del==0
  tau=[0 par(sys_tau)];
  m=length(tau)-1;
  xx=x;
  for j=1:m,
    xx=[xx x];  
  end;
else
  m=sys_ntau;
  xx=x;
  for j=1:m,
    t_tau(j)=sys_tau(j,xx,par);
    xx=[xx x];
  end;
  tau=[0 t_tau];
end;


for j=0:m
  A(:,:,j+1)=sys_deri(xx,par,j,[],[]);
end;

n1=zeros(size(l0));

for j=1:length(l0)

  corr=1+epsi;

  l=l0(j);

  D=l*eye(n);
  dD=eye(n);

  for jj=0:m
    B=A(:,:,jj+1)*exp(-l*tau(jj+1));
    D=D-B;
    dD=dD+tau(jj+1)*B;
  end;

  [E1,E2]=eig(D);

  [dummy,k]=min(abs(diag(E2)));

   v=E1(:,k);

  % choose normalisation c for v:

  vr=real(v);
  vi=imag(v);
  vn=vr'*vr+vi'*vi;

  c=v'/vn; % complex conjugate transpose

  while corr>epsi  

    D=l*eye(n);
    dD=eye(n);

    for jj=0:m
      B=A(:,:,jj+1)*exp(-l*tau(jj+1));
      D=D-B;
      dD=dD+tau(jj+1)*B;
    end;

    J=[D dD*v];
    J(n+1,:)=[c 0];

    r=-D*v;
    r(n+1,1)=-(c*v-1);

    d=J\r;

    l=l+d(n+1);
    v=v+d(1:n);

    n1(j)=n1(j)+1;

    corr=norm(d);

    if n1(j)>=max_n
      break;
    end;
   
  end;    

  if corr>epsi
    n1(j)=-1;
  end;

  l1(j)=l;

end;

return;

