function [J,res]=psol_jac(c,T,psol_prof,t,m,par,free_par,ph)

% function [J,res]=psol_jac(c,T,profile,t,deg,par,free_par,phase)
% INPUT:
%	c collocation parameters in [0,1]^deg
%	T period 
%	profile profile in R^(n x deg*l+1)
%	t representation points in [0,1]^(deg*l+1)
%	deg degree piecewise polynomial
%	par current parameter values in R^p
%	free_par free parameters numbers in N^d 
%	phase use phase condition or not (s = 1 or 0)
% OUTPUT: 
%	J jacobian in R^(n*deg*l+n+s x n*deg*l+1+n+d)
%	res residual in R^(n*deg*l+n+s)

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

n=size(psol_prof,1); % system dimension

tp_del=nargin('sys_tau');
if tp_del==0
  n_tau=sys_tau; % delay numbers
  tau=par(n_tau); % delay values
  tT=tau/T;
  d=length(n_tau); % number of delays
else
  d=sys_ntau; % number of delays
end;

l=(length(t)-1)/m; % number of intervals

% check:

if l~=floor(l)
  err=[length(t) m]
  error('PSOL_JAC: t does not contain l intervals of m points!');
end;

if length(c)~=m & ~isempty(c)
  err=[length(c) m]
  error('PSOL_JAC: wrong number of collocation parameters!');
end;

% init J, res:

mn=m*n;
nml=mn*l;
nml_n_1=nml+n+1;

J=zeros(nml+ph+n,nml_n_1+length(free_par));
res=zeros(nml+ph+n,1);

if isempty(c)
  c=poly_gau(m);
  gauss_c=c;
  non_gauss=0;
else
  gauss_c=poly_gau(m);
  non_gauss=1;
end;

% phase condition initialisation:

if ph
  
  gauss_abs=ones(1,m);
  g=poly_gau(m-1);
  for k=1:m
    for j=1:m-1
      gauss_abs(k)=gauss_abs(k)/(gauss_c(k)-g(j));
    end;
    for j=1:m
      if j~=k
        gauss_abs(k)=gauss_abs(k)/(gauss_c(k)-gauss_c(j));
      end;
    end;
  end;
  gauss_abs=gauss_abs/sum(gauss_abs);
end;

% more initialisation:

for m_i=1:m
  all_dPa(m_i,:)=poly_dla((0:m)/m,c(m_i));
  all_Pa(m_i,:)=poly_lgr((0:m)/m,c(m_i));
end;

% for all collocation points, make equation:

Pb=[];
x_tau=[];

for l_i=1:l

  % determine index for a in profile:

  index_a=(l_i-1)*m+1;

  i_index_a=(index_a-1)*n;

  t_start=t(index_a);
 
  hhh=t(index_a+m)-t_start;

  for m_i=1:m

    i=index_a+m_i-1;

    i_range=(i-1)*n+1:i*n;

    % determine c  

    col=t_start+c(m_i)*hhh;

    % determine dPa for a:

    dPa=all_dPa(m_i,:)/hhh;

    % add dPa for da in J:
    
    for k=0:m,
      kk=i_index_a+k*n;
      J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)+eye(n)*dPa(k+1);
    end;

    % add sum a*dPa to res:

    u_prime=psol_prof(:,index_a:index_a+m)*dPa';

    res(i_range,1)=u_prime; 

    %  determine Pa for a:

    Pa=all_Pa(m_i,:); 

    % phase_condition:
    
    if ph & ~non_gauss
      fup=gauss_abs(m_i)*hhh*u_prime';
      i_l_i=(l_i-1)*mn;
      for q=0:m
        qq=i_l_i+q*n;
        J(nml_n_1,qq+1:qq+n)=J(nml_n_1,qq+1:qq+n)+Pa(q+1)*fup;
      end;
    end;

    % compute x:

    x=psol_prof(:,index_a:index_a+m)*Pa';

    % compute tau, c_tau, x_tau, dx_tau
    xx=x;

    for tau_i=1:d

      if tp_del~=0
        tT(tau_i)=sys_tau(tau_i,xx,par)/T;
      end;

      c_tau_i=col-tT(tau_i);
      while c_tau_i<0, 
        c_tau_i=c_tau_i+1; 
      end; 
 
      c_tau(tau_i)=c_tau_i;      

      % determine index for b in profile:

      index_b_i=length(t)-m;
      while (c_tau_i<t(index_b_i)) 
        index_b_i=index_b_i-m; 
      end;
   
      index_b(tau_i)=index_b_i;

      % c transformed to [0,1] and hhh_tau
     
      hhh_tau(tau_i)=t(index_b_i+m)-t(index_b_i);
      c_tau_trans(tau_i)=(c_tau_i-t(index_b_i))/hhh_tau(tau_i);

      % determine Pb for b:

      Pb(tau_i,:)=poly_elg(m,c_tau_trans(tau_i));

      % compute x_tau:
  
      x_tau(:,tau_i)=psol_prof(:,index_b_i:index_b_i+m)*Pb(tau_i,:)';

      dPb(tau_i,:)=poly_del(m,c_tau_trans(tau_i))/hhh_tau(tau_i);
      dx_tau(:,tau_i)=psol_prof(:,index_b_i:index_b_i+m)*dPb(tau_i,:)';
      xx=[x x_tau];

    end;

    % determine A0:
  
    T_A0=T*sys_deri(xx,par,0,[],[]);

    % add -T*A0*Pa for da in J:

    for k=0:m
      kk=i_index_a+k*n; 
      J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)-T_A0*Pa(k+1);  
    end;

    % determine parameter derivatives:

    for p_i=1:length(free_par)
      df=sys_deri(xx,par,[],free_par(p_i),[]);
      J(i_range,nml_n_1+p_i)=-T*df; 
    end;  

    % compute f(x,x_tau):

    f=sys_rhs(xx,par);

    % add -f for dT in J:

    J(i_range,nml_n_1)=J(i_range,nml_n_1)-f;

    % add Tf in res:

    res(i_range,1)=res(i_range,1)-T*f;

    TPb=T*Pb;

    for t_i=1:d

      % determine A1:

      A1=sys_deri(xx,par,t_i,[],[]);

      % add -T*A1*Pb for db and A1*dx_tau*dtau*Pa for da in J:

      i_index_b=(index_b(t_i)-1)*n;
      for k=0:m
        kk=i_index_b+n*k;
        J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)-A1*TPb(t_i,k+1);  
      end;

      % add -T*A1*sum b*dP*dc_tau for dT in J:

      J(i_range,nml_n_1)=J(i_range,nml_n_1)-A1*dx_tau(:,t_i)*tT(t_i);

      if tp_del~=0
        dtau=sys_dtau(t_i,xx,par,0,[]);
        % add A1*dx_tau*dtau*Pa for da in J 
        for k=0:m
          kk=i_index_a+k*n;
          J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n) ...
                               +A1*(dx_tau(:,t_i)*dtau)*Pa(k+1);
        end;

        % delay function depends on delayed terms
        % (e.g. tau_2(...)=g(x(t),x(t-tau_1)), tau_1 - constant )

        for t_ii=1:t_i-1
          dtau=sys_dtau(t_i,xx,par,t_ii,[]);
          A1_dx_dtau=A1*(dx_tau(:,t_i)*dtau);
          i_index_b=(index_b(t_ii)-1)*n;
        % add A1*dx_tau*dtau*Pb for db in J
          for k=0:m
            kk=i_index_b+n*k;
            J(i_range,kk+1:kk+n)=J(i_range,kk+1:kk+n)+A1_dx_dtau*Pb(t_ii,k+1);
          end;
        % for dT in J
            J(i_range,nml_n_1)=J(i_range,nml_n_1) ...
                              +A1_dx_dtau*dx_tau(:,t_ii)*tT(t_ii)/T;
        % tau_1 - continuation parameter
           if norm(dtau)~=0 
             for p_i=1:length(free_par)
               dtau_fp=sys_dtau(t_ii,xx,par,[],free_par(p_i));
               if dtau_fp~=0 
                 J(i_range,nml_n_1+p_i)=J(i_range,nml_n_1+p_i) ...
                                -A1_dx_dtau*(dx_tau(:,t_ii))/T; 
               end;
             end;
           end;
        end;
      end;

      for p_i=1:length(free_par)
        if tp_del==0 & free_par(p_i)==n_tau(t_i)  
          J(i_range,nml_n_1+p_i)=J(i_range,nml_n_1+p_i)+A1*dx_tau(:,t_i);
        elseif tp_del~=0
          dtau=sys_dtau(t_i,xx,par,[],free_par(p_i));
          J(i_range,nml_n_1+p_i)=J(i_range,nml_n_1+p_i) ...
              +A1*(dx_tau(:,t_i)*dtau);
        end;
      end;
    end;
  end;
end;

% periodicity condition:

J(nml+1:nml+n,1:n)=eye(n);
J(nml+1:nml+n,nml+1:nml+n)=-eye(n);
res(nml+1:nml+n,1)=psol_prof(:,1)-psol_prof(:,size(psol_prof,2));

% phase condition:

if ph & non_gauss,
  for l_i=1:l
    index_a=(l_i-1)*m+1;
    for k=1:m
      fac=gauss_abs(k)*(t((l_i-1)*m+1)-t(l_i*m+1));
      dPa=poly_dla(t(index_a:index_a+m),gauss_c(k));
      Pa=poly_lgr(t(index_a:index_a+m),gauss_c(k));
      u_prime=psol_prof(:,index_a:index_a+m)*dPa';
      for q=1:m+1
        J(nml_n_1,(l_i-1)*mn+1+(q-1)*n:(l_i-1)*mn+q*n)= ...
	  J(nml_n_1,(l_i-1)*mn+1+(q-1)*n:(l_i-1)*mn+q*n) + fac*Pa(q)*u_prime';
      end;
    end;
  end;
end;

if ph
  res(nml_n_1,1)=0;
end;

return;
