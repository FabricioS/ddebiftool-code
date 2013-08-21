function [eig_mesh,eig_v]=mult_dbl(period,profile,mesh,degree,col,par)

% function eig_v=mult_dbl(period,profile,mesh,degree,col,par)
% INPUT: 
%	period period of solution
%	profile periodic solution profile
%       mesh periodic solution mesh (if empty, mesh is assumed uniform)
%	degree degree of piecewise polynomials
%	col collocation points
%       par current parameter values in R^p
% OUTPUT:
%	eig_mesh mesh for eigenvector
%       eig_v eigenvector corresponding to multiplier near -1

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

% is mesh is empty assume equidistant mesh:
if isempty(mesh)
  mesh=0:1/(size(profile,2)-1):1;
end;

m=degree;

tp_del=nargin('sys_tau');
if tp_del==0
  tau=par(sys_tau);
  tT=max(tau)/period;
  d=length(tau);
else
  d=sys_ntau;
  tot_tau=[];
  for t_m=1:length(mesh)
    x=profile(:,t_m);
    xx=x;
    for j=1:d
      tau=sys_tau(j,xx,par);
      tot_tau=[tot_tau tau];
      tau=tau/period;
      t_tau=mesh(t_m)-tau;
      while t_tau<0,
        t_tau=t_tau+1;
      end;
      index_b=length(mesh)-m;
      while (t_tau<mesh(index_b))
        index_b=index_b-m;
      end;
      hhh_tau=mesh(index_b+m)-mesh(index_b);
      t_tau_trans=(t_tau-mesh(index_b))/hhh_tau;
      Pb=poly_elg(m,t_tau_trans);
      x_tau=profile(:,index_b:index_b+m)*Pb';
      xx=[xx x_tau];
    end;
  end;
  tT=max(tot_tau)/period;
end;

k=floor(tT); 

for l=length(mesh):-m:1
  if mesh(l)<k+1-tT
    break;
  end;
end;

ext_mesh=mesh(l:length(mesh))-k-1;

for l=k:-1:0
  ll=length(ext_mesh);
  ext_mesh(ll:ll-1+length(mesh))=mesh-l;
end;

T=period;
t=mesh;
tt=ext_mesh;
if isempty(col)
  co=poly_gau(degree);
else 
  co=col;
end;

n=size(profile,1);
l=(length(t)-1)/m;
mm=length(co);
ll=(length(tt)-1)/mm;

if l~=floor(l)
  err=[length(t) m]
  error('MULT_DBL: t does not contain l intervals of m points!');
end;
if ll~=floor(ll)
  err=[length(tt) mm]
  error('MULT_DBL: tt does not contain ll intervals of m points!');
end;

% init J:

k=1;
while tt(k)<0,
  k=k+1;
end;
J=zeros(n*(length(tt)-k),n*(mm*ll+1));

% more initialisation:

for m_i=1:mm
  all_dPa(m_i,:)=poly_dla((0:mm)/mm,co(m_i));
  all_Pa(m_i,:)=poly_lgr((0:mm)/mm,co(m_i));
end;

% for all collocation points, make equation:

i=0;

for ll_i=1:ll

  if tt((ll_i-1)*mm+1)>=0
    
    % determine index for a in profile:

    index_a=(ll_i-1)*mm+1;

    hhh=tt(index_a+mm)-tt(index_a);

    for m_i=1:mm

      i=i+1;

      i_range=(i-1)*n+1:i*n;

      cc=tt(index_a)+co(m_i)*hhh;
     
      % determine dPa for a:

      dPa=all_dPa(m_i,:)/hhh;

      % add dPa for da in J:

      for k=1:mm+1,
        kk=(index_a-1)*n+(k-1)*n;
        J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)+eye(n)*dPa(k); 
      end;

      %  determine Pa for a:

      Pa=all_Pa(m_i,:); 

      % compute x:

      x=psol_eva(profile,t,cc,m);

      % compute tau, c_tau, x_tau, dx_tau
      xx=x;

      for t_i=1:d

      if tp_del~=0
        tau(t_i)=sys_tau(t_i,xx,par);
      end;

        % determine c_tau:

        c_tau(t_i)=cc-tau(t_i)/T;

        % determine index for b in profile:

        index_b(t_i)=length(tt)-mm;
        while (c_tau(t_i)<tt(index_b(t_i))) 
          index_b(t_i)=index_b(t_i)-mm;
        end;

        % c transformed to [0,1] and hhh_tau
     
        hhh_tau=tt(index_b(t_i)+mm)-tt(index_b(t_i));
        c_tau_trans=(c_tau(t_i)-tt(index_b(t_i)))/hhh_tau;

        % determine Pb for b:

        Pb(t_i,:)=poly_elg(m,c_tau_trans);

        % compute x_tau:

        x_tau(:,t_i)=psol_eva(profile,t,c_tau(t_i),m);

        xx=[x x_tau];

      end;

      % determine A0:

      A0=sys_deri(xx,par,0,[],[]);

      % add -T*A0*Pa for da in J:

      for k=1:m+1
         kk=(index_a-1)*n+(k-1)*n;
         J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)-T*A0*Pa(k);
      end;

      if tp_del~=0
        for j=1:d
          xxx=x_tau(:,j);
          for jj=1:d
            tau_sh=sys_tau(jj,xxx,par);
            c_tau_sh=c_tau(j)-tau_sh/T;
            x_tau_sh=psol_eva(profile,t,c_tau_sh,m);
            xxx=[xxx x_tau_sh];
          end;
          dx_tau(:,j)=T*sys_rhs(xxx,par);
          xxx=[];
        end;
      end;

      for t_i=1:d

        % determine A1:

        A1=sys_deri(xx,par,t_i,[],[]);

        % add -T*A1*Pb for db in J:

        for k=1:m+1
          kk=(index_b(t_i)-1)*n+(k-1)*n;
          J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)-T*A1*Pb(t_i,k);
        end;

        if tp_del~=0
          dtau=sys_dtau(t_i,xx,par,0,[]);
          % add A1*dx_tau*dtau*Pa for da in J:
          for k=1:m+1
            kk=(index_a-1)*n+(k-1)*n;
            J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n) ...
                                 +A1*(dx_tau(:,t_i)*dtau)*Pa(k); 
          end;

          % delay function depends on delayed terms
          % (e.g. tau_2(...)=g(x(t),x(t-tau_1)), tau_1 - constant )

          for t_ii=1:t_i-1
            dtau=sys_dtau(t_i,xx,par,t_ii,[]);
            A1_dx_dtau=A1*(dx_tau(:,t_i)*dtau);
            % add A1*dx_tau*dtau*Pb for db in J
            for k=1:m+1
              kk=(index_b(t_ii)-1)*n+(k-1)*n;
              J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)+A1_dx_dtau*Pb(t_ii,k);
            end;
          end;
        end;

      end;

    end;

  end;

end;

% lower triangularize:

k=1;
while tt(k)<0,
  k=k+1;
end;

for i=1:(size(J,2)-k*n)/(mm*n)
  X=J(mm*n*(i-1)+1:mm*n*i,k*n+1+mm*n*(i-1):k*n+mm*n*i);
  [L,U,P]=lu(X);
  J(mm*n*(i-1)+1:mm*n*i,:)=inv(L)*P*J(mm*n*(i-1)+1:mm*n*i,:);
  for j=2:mm*n
    for q=1:j-1
      J(mm*n*(i-1)+j,k*n+q+mm*n*(i-1))=0;
    end;
  end;
end;

for i=1:(size(J,2)-k*n)/(mm*n),
  for j=1:mm*n,
    pivot=J(mm*n*(i-1)+j,k*n+mm*n*(i-1)+j);
    for q=mm*n*(i-1)+j+1:size(J,1),
      if J(q,k*n+mm*n*(i-1)+j)~=0,
        J(q,:)=J(q,:)-J(q,k*n+mm*n*(i-1)+j)*J(mm*n*(i-1)+j,:)/pivot;
        J(q,k*n+mm*n*(i-1)+j)=0;
      end;
    end;
  end;
end;

% extract M:

if k*n<size(J,2)-k*n+1
  M0=J(size(J,1)-k*n+1:size(J,1),size(J,2)-k*n+1:size(J,2));
  M1=J(size(J,1)-k*n+1:size(J,1),1:k*n);
  M=-inv(M0)*M1;
else
  overlap=k*n-(size(J,2)-k*n+1)+1;
  M0=J(size(J,1)-k*n+1+overlap:size(J,1),size(J,2)-k*n+1+overlap:size(J,2));
  M1=J(size(J,1)-k*n+1+overlap:size(J,1),1:k*n);
  M(1:overlap,k*n-overlap+1:k*n)=eye(overlap);
  M(overlap+1:k*n,1:k*n)=-inv(M0)*M1;
end;

[S1,S2]=eig(M);

s=diag(S2);

for i=1:length(s)
  if abs(imag(s(i)))>0
    s(i)=0;
  end;
end;

[i1,i2]=min(abs(s+1));

if s(i2)==0
  err=s(i2)
  error('MULT_DBL: Could not find a candidate -1 multiplier.');
end;

V=S1(:,i2);

l1=size(J,1);
l0=l1/n;
l2=size(J,2);
l3=length(V);

vv=J(1:l1,l3+1:l2)\(-J(1:l1,1:l3)*V);

eig_v(1:n,1)=-vv(l1-n+(1:n));

for j=1:l0
  eig_v(1:n,j+1)=vv((j-1)*n+(1:n));
end;

for j=1:l0
  eig_v(1:n,l0+j+1)=-vv((j-1)*n+(1:n));
end;

eig_mesh=[t/2 0.5*(1+t(2:l0+1))];

return;
