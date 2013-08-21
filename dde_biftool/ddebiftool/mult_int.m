function M=mult_int(T,profile,t,m,tt,co,par,d_ac)

% function M=mult_int(T,profile,t,m,tt,co,par,d_ac)
% INPUT:
%	T period
%	profile in R^(n x m*l+1)
%	t representation points in [0,1]^(m*l+1)
%	m order polynomial
%	tt extended stability time mesh, t in tt
%	co nonempty row of stability collocation parameters
%	par current parameter values
%       d_ac (only for state-dependent delays) tau<d_ac is treated as 
%             tau<0 (stability is not computed)
% OUTPUT: 
%	M approximation of monodromy matrix

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

n=size(profile,1);
l=(length(t)-1)/m;
mm=length(co);
ll=(length(tt)-1)/mm;

lttmm=length(tt)-mm;

tp_del=nargin('sys_tau');
if tp_del==0
  tau=par(sys_tau);
  d=length(tau);
else
  d=sys_ntau;
end;

if l~=floor(l)
  err=[length(t) m]
  error('MULT_INT: t does not contain l intervals of m points!');
end;
if ll~=floor(ll)
  err=[length(tt) mm]
  error('MULT_INT: tt does not contain ll intervals of m points!');
end;

% init J:

mmn=mm*n;

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

Pb=[];
x_tau=[];

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
        if (tau(t_i)<d_ac),
          M=[];
          s=strcat('WARNING: delay number_',num2str(t_i),' is negative, no stability computed.');
          disp(s);
          return;
        end;
      end;

        % determine c_tau:

        c_tau=cc-tau(t_i)/T;
        c_tau_i(t_i)=c_tau;

        % determine index for b in profile:

        index_b_i=lttmm;
        while (c_tau<tt(index_b_i)) 
          index_b_i=index_b_i-mm;
        end;

        index_b(t_i)=index_b_i;

        % c transformed to [0,1] and hhh_tau
     
        hhh_tau=tt(index_b_i+mm)-tt(index_b_i);
        c_tau_trans=(c_tau-tt(index_b_i))/hhh_tau;

        % determine Pb for b:

        Pb(t_i,:)=poly_elg(m,c_tau_trans);

        % compute x_tau:

        x_tau(:,t_i)=psol_eva(profile,t,c_tau,m);

        xx=[x x_tau];

      end;

      % determine A0:

      A0=sys_deri(xx,par,0,[],[]);

      % add -T*A0*Pa for da in J:

      for k=1:m+1
        kk=(index_a-1)*n+(k-1)*n;
        J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)-T*A0*Pa(k);  
      end;

      TPb=T*Pb;

      if tp_del~=0
        % compute time derivative of solution 
        for j=1:d
          xxx=x_tau(:,j);
          for jj=1:d
            tau_sh=sys_tau(jj,xxx,par);
            c_tau_sh=c_tau_i(j)-tau_sh/T;
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
          J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)-A1*TPb(t_i,k);  
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

kn=k*n;

for i=1:(size(J,2)-kn)/mmn
  X=J(mmn*(i-1)+1:mmn*i,kn+1+mmn*(i-1):kn+mmn*i);
  [L,U,P]=lu(X);
  J(mmn*(i-1)+1:mmn*i,:)=L\P*J(mmn*(i-1)+1:mmn*i,:);
  for j=2:mmn
    for q=1:j-1
      J(mmn*(i-1)+j,kn+q+mmn*(i-1))=0;
    end;
  end;
end;

lJ1=size(J,1);
lJ2=size(J,2);

for i=1:(lJ2-kn)/mmn,
  mmn_i_1=mmn*(i-1);
  kn_mmn_i_1=kn+mmn*(i-1);
  for j=1:mmn,
    kn_mmn_i_1_j=kn_mmn_i_1+j;
    pivot=J(mmn_i_1+j,:)/J(mmn_i_1+j,kn_mmn_i_1_j);
    qq=mmn_i_1+j+1:lJ1;
    ql=qq(J(qq,kn_mmn_i_1_j)~=0);
    for l=1:length(ql)
      q=ql(l);
      J(q,:)=J(q,:)-J(q,kn_mmn_i_1_j)*pivot;
      J(q,kn_mmn_i_1_j)=0;
    end;
  end;
end;

% extract M:

if kn<size(J,2)-kn+1
  M0=J(lJ1-kn+1:lJ1,lJ2-kn+1:lJ2);
  M1=J(lJ1-kn+1:lJ1,1:kn);
  M=-M0\M1;
else
  overlap=kn-(lJ2-kn+1)+1;
  M0=J(lJ1-kn+1+overlap:lJ1,lJ2-kn+1+overlap:lJ2);
  M1=J(lJ1-kn+1+overlap:lJ1,1:kn);
  M(1:overlap,kn-overlap+1:kn)=eye(overlap);
  M(overlap+1:kn,1:kn)=-M0\M1;
end;

return;
