function [point,success]=p_correc(point,free_par,step_cnd,method,p_nr,...
                                  previous,d_nr,tz)

% function [point,success]=p_correc(point0,free_par,step_cnd,method,adapt,
%                                 previous)
% INPUT:
%   point0 initial point guess
%   free_par free parameter numbers in N^d 
%       step_cnd steplength condition(s) as point(s)
%   method method parameters
%   adapt if zero or absent, do not adapt mesh; if one, always adapt 
%       previous (optional) previously computed branch point (used, in case 
%            of periodic solutions or connecting orbits, to
%            minimize phase shift)
% OUTPUT:
%       point final solution
%       success nonzero for successfull computation
%
% OPTIONAL EXTRA INPUT:
%       d_nr number of delay crossing zero
%       tz (for psol only) time point for which tau(tz)=0 and dtau(tz)/dt=0

% (c) DDE-BIFTOOL v. 2.02, 16/6/2002

% optional parameters:

if ~exist('p_nr'),
  p_nr=0;
end;
if ~exist('previous'),
  previous=[];
end;

if ~exist('d_nr'),
  d_nr=0;
else
  d_tau=sys_ntau;
end;

% some method parameters:

max_iter=method.newton_max_iterations; % max number of newton iterations
nmon_iter=method.newton_nmon_iterations; % max number of nonmonotone iterations
conv_r=method.halting_accuracy; % required accuracy for convergence
print_r=method.print_residual_info; % print residual evolution

% remove old stability info if present:

if isfield(point,'stability')
  point.stability=[];
end;

% initialize:

switch point.kind,
  case 'stst',
    n=size(point.x,1);
    p_start=n;
  case 'fold',
    n=size(point.x,1);
    p_start=2*n;
    c=point.v'/(point.v'*point.v);
  case 'hopf',
    n=size(point.x,1);
    p_start=3*n+1;
    vr=real(point.v);
    vi=imag(point.v);
    vn=vr'*vr+vi'*vi;
    c=point.v'/vn; % complex conjugate transpose
  case {'psol','hcli'},
    n=size(point.profile,1);
    col=method.collocation_parameters;
    mm=point.degree;
    if ~isempty(col) & length(col)~=mm
      err=[point.degree mm]
      error('P_CORREC: number of collocation parameters differs from degree!');
    end;
    ll=floor((size(point.profile,2)-1)/mm);
    if ll~=(size(point.profile,2)-1)/mm,
      error('P_CORREC: point mesh does not contain l intervals of m points!');
    end;
    p_start=n*mm*ll+1+n;
    if d_nr==0
      ph=method.phase_condition;
    else
      ph=0;
    end;
    if isempty(point.mesh)
      mesh=0:1/(ll*mm):1;
    else
      mesh=point.mesh;
      ma=method.adapt_mesh_before_correct;
      if p_nr>1 & mod(p_nr,ma)==0 
        % do not adapt mesh when p_nr=0
        % do not (yet) adapt mesh when p_nr=1
        % adapt mesh when p_nr>1 & mod(p_nr,ma)=0
        new_mesh=psol_msh(mesh,mm,point.profile,ll,point.degree);
        switch point.kind,
          case 'psol',
            point.profile=psol_eva(point.profile,mesh,new_mesh,mm);
          case 'hcli',
            point.profile=hcli_eva(point.profile,mesh,new_mesh,mm);
        end;
        point.mesh=new_mesh;
        mesh=new_mesh;
      end;
    end;
  otherwise,
    err=point.kind
    error('P_CORREC: point kind not recognized.');
  end;

% store parameters:

par=point.parameter;

% perform Newton-Raphson iterations:

for i=1:max_iter

  % compute & apply corrections:

  switch point.kind  

    case 'stst',
      % produce jacobian
      [J,res]=stst_jac(point.x,par,free_par);
      % add (linear) steplength conditions
      for j=1:length(step_cnd)
        J(n+j,1:n)=step_cnd(j).x';
        for k=1:length(free_par)
          J(n+j,n+k)=step_cnd(j).parameter(free_par(k));
        end;
        res(n+j,1)=0;
      end;
      if d_nr~=0
        % add condition on zero delay
        kk=size(J,1);
        for jj=0:d_tau
          xx(:,jj+1)=point.x;
        end;
        J(kk+1,1:n)=0;
        for j=0:d_tau
          J(kk+1,1:n)=J(kk+1,1:n)+sys_dtau(d_nr,xx,par,j,[]);
        end;
        for k=1:length(free_par)
          J(kk+1,n+k)=sys_dtau(d_nr,xx,par,[],free_par(k));
        end;
        res(kk+1,1)=sys_tau(d_nr,xx,par);
      end;
      % add extra conditions
      if method.extra_condition
        [resi,condi]=sys_cond(point);
        for j=1:length(condi)
          k=size(J,1);
          res(k+1)=resi(j);
          J(k+1,1:n)=condi(j).x';
          for o=1:length(free_par)
            J(k+1,n+o)=condi(j).parameter(free_par(o));
          end;
        end;
      end;
      % solve linear system
      if size(J,1)~=size(J,2)
        disp('P_CORREC warning: use of nonsquare Jacobian.');
      end;

      dx=J\res;
      % apply non-parameter corrections
      point.x=point.x-dx(1:n);

    case 'fold',
      % produce jacobian
      [J,res]=fold_jac(point.x,point.v,par,free_par,c);
      % add (linear) steplength conditions
      for j=1:length(step_cnd)
        J(2*n+1+j,1:n)=step_cnd(j).x';
        J(2*n+1+j,n+1:2*n)=step_cnd(j).v';
        for k=1:length(free_par)
          J(2*n+1+j,2*n+k)=step_cnd(j).parameter(free_par(k));
        end;
        res(2*n+1+j,1)=0;
      end;
      if d_nr~=0
        % add condition on zero delay
        kk=size(J,1);
        for jj=0:d_tau
          xx(:,jj+1)=point.x;
        end;
        J(kk+1,1:2*n+2)=0;
        for j=0:d_tau
          J(kk+1,1:n)=J(kk+1,1:n)+sys_dtau(d_nr,xx,par,j,[]);
        end;
        for k=1:length(free_par)
          J(kk+1,2*n+k)=sys_dtau(d_nr,xx,par,[],free_par(k));
        end;
        res(kk+1,1)=sys_tau(d_nr,xx,par);
      end;
      % add extra conditions
      if method.extra_condition
        [resi,condi]=sys_cond(point);
        for j=1:length(condi)
          k=size(J,1);
          res(k+1)=resi(j);
          J(k+1,1:n)=condi(j).x';
          J(k+1,n+1:2*n)=condi(j).v';
          for o=1:length(free_par)
            J(k+1,2*n+o)=condi(j).parameter(free_par(o));
          end;
        end; 
      end;
      % solve linear system
      if size(J,1)~=size(J,2)
        disp('P_CORREC warning: use of nonsquare Jacobian.');
      end;
      dx=J\res;
      % apply non-parameter corrections
      point.x=point.x-dx(1:n);
      point.v=point.v-dx(n+1:2*n);

    case 'hopf'  
      % produce jacobian
      [J,res]=hopf_jac(point.x,point.omega,point.v,par,free_par,c);
        % add (linear) steplength conditions
        for j=1:length(step_cnd)
          J(3*n+2+j,1:n)=step_cnd(j).x';
          J(3*n+2+j,n+1:2*n)=real(step_cnd(j).v)';
          J(3*n+2+j,2*n+1:3*n)=imag(step_cnd(j).v)';
          J(3*n+2+j,3*n+1)=step_cnd(j).omega;
          for k=1:length(free_par)
            J(3*n+2+j,3*n+1+k)=step_cnd(j).parameter(free_par(k));
          end;
          res(3*n+2+j)=0;
        end;
      if d_nr~=0
        % add condition on zero delay
        kk=size(J,1);
        for jj=0:d_tau
          xx(:,jj+1)=point.x;
        end;
        J(kk+1,1:3*n+3)=0;
        for j=0:d_tau
          J(kk+1,1:n)=J(kk+1,1:n)+sys_dtau(d_nr,xx,par,j,[]);
        end;
        for k=1:length(free_par)
          J(kk+1,3*n+1+k)=sys_dtau(d_nr,xx,par,[],free_par(k));
        end;
        res(kk+1,1)=sys_tau(d_nr,xx,par);
      end;
      % add extra conditions
      if method.extra_condition
        [resi,condi]=sys_cond(point);
        for j=1:length(condi)
          k=size(J,1);
          res(k+1)=resi(j);
          J(k+1,1:n)=condi(j).x';
          J(k+1,n+1:2*n)=real(condi(j).v');
          J(k+1,2*n+1:3*n)=imag(condi(j).v');
          J(k+1,3*n+1)=condi(j).omega;
          for o=1:length(free_par)
            J(k+1,3*n+1+o)=condi(j).parameter(free_par(o));
          end;
        end;
      end;
      % solve linear system
      if size(J,1)~=size(J,2)
        disp('P_CORREC warning: use of nonsquare Jacobian.');
      end;
      dx=J\res;
      % apply non-parameter corrections
      point.x=point.x-dx(1:n);
      point.v=point.v-dx(n+1:2*n)-sqrt(-1)*dx(2*n+1:3*n);
      point.omega=point.omega-dx(3*n+1);

    case 'psol',
      % produce jacobian
      [J,res]=psol_jac(col,point.period,point.profile,mesh,mm,par,...
                         free_par,ph);
      % add (linear) steplength conditions
      for j=1:length(step_cnd)
        for k=1:ll*mm+1
          J(n*(mm*ll+1)+ph+j,(k-1)*n+1:k*n)=step_cnd(j).profile(:,k)';
        end;
        J(n*(mm*ll+1)+ph+j,n*(mm*ll+1)+1)=step_cnd(j).period;
        for k=1:length(free_par)
          J(n*(mm*ll+1)+ph+j,n*(mm*ll+1)+1+k)=...
                                  step_cnd(j).parameter(free_par(k));
        end;
        res(n*(mm*ll+1)+ph+j)=0;
      end;
      if d_nr~=0
        %%%% add 2 conditions: on zero delay (1) and 
        %%%% on zero time derivative of delay (2)

        T=point.period;

        % compute x_tau, dx_tau, d2x_tau 

        for tau_i=1:d_tau+1
          if tau_i==1
            tau(tau_i)=0;
          else
            tau(tau_i)=sys_tau(tau_i-1,x_tau,par)/T;
          end;
          tz_tau=tz-tau(tau_i);
          while tz_tau<0
            tz_tau=tz_tau+1;
          end;

          % determine index for b in profile:
          index_b_i=length(mesh)-mm;
          while (tz_tau<mesh(index_b_i))
            index_b_i=index_b_i-mm; 
          end;
          index_b(tau_i)=index_b_i;

          % tz_tau transformed to [0,1] and hhh
          hhh=mesh(index_b_i+mm)-mesh(index_b_i);
          tz_tau_trans=(tz_tau-mesh(index_b_i))/hhh;

          Pb(tau_i,:)=poly_elg(mm,tz_tau_trans);
          x_tau(:,tau_i)=point.profile(:,index_b_i:index_b_i+mm)*Pb(tau_i,:)';
          dPb(tau_i,:)=poly_del(mm,tz_tau_trans)/hhh;
          dx_tau(:,tau_i)=point.profile(:,index_b_i:index_b_i+mm)*...
                                                    dPb(tau_i,:)';
          d2Pb(tau_i,:)=poly_d2l(mm,tz_tau_trans)/(hhh^2);
          d2x_tau(:,tau_i)=point.profile(:,index_b_i:index_b_i+mm)*...
                                                    d2Pb(tau_i,:)';
        end;

        % fill in J and res 

        kj=size(J,1);    
        rr=n*(mm*ll+1);
        J(kj+1:kj+2,1:rr+1+length(free_par))=0;
        res(kj+1:kj+2,1)=0;
    
        for tau_i=1:d_tau+1
          dtau=sys_dtau(d_nr,x_tau,par,tau_i-1,[]);
          sd=0;
          for t_i=1:d_tau+1
            d2tau=sys_dtau(d_nr,x_tau,par,[tau_i-1 t_i-1],[])*dx_tau(:,t_i);
            sd=sd+d2tau;
            % cond. (2), add for dT in J
            J(kj+2,rr+1)=J(kj+2,rr+1)+ ...
                         (d2tau'*dx_tau(:,tau_i))*tau(t_i)/T;
          end; 
          i_index_b=(index_b(tau_i)-1)*n;
          for k=0:mm
            kk=i_index_b+n*k;
            % cond. (1), add dtau/dx, dtau/dx_tau, ...
            J(kj+1,kk+1:kk+n)=J(kj+1,kk+1:kk+n)+dtau*Pb(tau_i,k+1); 
            % cond. (2), add for time derivative of delay
            J(kj+2,kk+1:kk+n)=J(kj+2,kk+1:kk+n)+dtau*dPb(tau_i,k+1)+ ...
                                                (sd)'*Pb(tau_i,k+1);
          end;

          % cond. (1), add for dT in J
          J(kj+1,rr+1)=J(kj+1,rr+1)+dtau*dx_tau(:,tau_i)*tau(tau_i)/T;

          % cond. (2), add for dT in J
          J(kj+2,rr+1)=J(kj+2,rr+1)+dtau*d2x_tau(:,tau_i)*tau(tau_i)/T;      

          % cond. (2), derivatives wrt free parameter
          for k=1:length(free_par)
            J(kj+2,rr+1+k)=J(kj+2,rr+1+k)+ ...
             sys_dtau(d_nr,x_tau,par,tau_i-1,free_par(k))*dx_tau(:,tau_i);
          end;

          % conds. (1) and (2), if constant delay is free parameter

          if tau_i>1 & norm(dtau)~=0, %delay d_nr depends on x(t-tau(tau_i-1))
            for k=1:length(free_par)
              dtau_fp=sys_dtau(tau_i-1,x_tau,par,[],free_par(k));
              if dtau_fp~=0 % = tau(tau_i-1) is free parameter
                J(kj+1,rr+1+k)=J(kj+1,rr+1+k)-dtau*dx_tau(:,tau_i)/T;
                for t_i=1:d_tau+1
                  J(kj+2,rr+1+k)=J(kj+2,rr+1+k)- ...             
                    (sys_dtau(d_nr,x_tau,par,[t_i-1 tau_i-1],[])* ...
                               dx_tau(:,tau_i))'*(dx_tau(:,t_i))/T;
                  if t_i==tau_i
                    J(kj+2,rr+1+k)=J(kj+2,rr+1+k)-dtau*d2x_tau(:,tau_i)/T;
                  end;
                end;
              end;
            end;
          end;

          % cond. (2), residual 
          res(kj+2,1)=res(kj+2,1)+dtau*dx_tau(:,tau_i);
        end;

        % cond. (1), derivatives wrt free paramter
        for k=1:length(free_par)
          J(kj+1,rr+1+k)=J(kj+1,rr+1+k)+...
                                sys_dtau(d_nr,x_tau,par,[],free_par(k));
        end;

        % cond. (1), residual 
        res(kj+1,1)=sys_tau(d_nr,x_tau,par);

        if norm(J(kj+2,:))==0 %dtau/dt=0 for constant tau
          ph=method.phase_condition;
          [J,res]=psol_jac(col,point.period,point.profile,mesh,mm,par,...
                         free_par,ph); 
          J(kj+2,1:rr+1+length(free_par))=0;
          for k=1:length(free_par)
            J(kj+2,rr+1+k)=sys_dtau(d_nr,x_tau,par,[],free_par(k));
          end;
          res(kj+2,1)=sys_tau(d_nr,x_tau,par); 
        end;
      end;     
      % add extra conditions
      if method.extra_condition
        [resi,condi]=sys_cond(point);
        for j=1:length(condi)
          k=size(J,1);
          res(k+1)=resi(j);
          for o=1:ll*mm+1
            J(k+1,(o-1)*n+1:o*n)=condi(j).profile(:,o)';
          end;
          J(k+1,n*mm*ll+n+1)=condi(j).period;
          for o=1:length(free_par)
            J(k+1,n*mm*ll+n+1+o)=condi(j).parameter(free_par(o));
          end;
        end; 
      end;
      % solve linear system
      if size(J,1)~=size(J,2)
        disp('P_CORREC warning: use of nonsquare Jacobian.');
      end;
      dx=J\res;
      % apply non-parameter corrections
      for k=1:ll*mm+1
        point.profile(:,k)=point.profile(:,k)-dx((k-1)*n+1:k*n);
      end;
      point.period=point.period-dx(n*mm*ll+1+n);

    case 'hcli',
      % produce jacobian
      [J,res]=hcli_jac(col,point.period,point.profile,...
                       mesh,mm,par,free_par,ph,...
                       point.lambda_v,point.lambda_w,point.v,...
                       point.w,point.alpha,point.epsilon,point.x1,...
                       point.x2,previous);           
      % add (linear) steplength conditions
      for j=1:length(step_cnd)
        sJ=size(J,1);      
        for k=1:ll*mm+1
          J(sJ+1,(k-1)*n+1:k*n)=step_cnd(j).profile(:,k)';
        end;       
        J(sJ+1,n*(mm*ll+1)+1)=step_cnd(j).period;
        nb_par=length(free_par); 
        for k=1:nb_par
          J(sJ+1,n*(mm*ll+1)+1+k)=step_cnd(j).parameter(free_par(k));
        end;      
        J(sJ+1,n*(mm*ll+1)+1+nb_par+(1:n))=step_cnd(j).x1';
        J(sJ+1,n*(mm*ll+1)+nb_par+n+1+(1:n))=step_cnd(j).x2';      
        s1=length(step_cnd(j).v(1,:));
        for k=1:s1
          J(sJ+1,n*(mm*ll+1)+1+nb_par+2*n+(k-1)*n+(1:n))=...
              (step_cnd(j).v(:,k))';
        end;      
        J(sJ+1,n*(mm*ll+1)+1+2*n+nb_par+s1*n+(1:s1))=step_cnd(j).lambda_v';

        s2=size(step_cnd(j).w,2);
        for k=1:s2
         J(sJ+1,n*(mm*ll+1)+1+nb_par+2*n+s1*(n+1)+(k-1)*n+(1:n))=...
              step_cnd(j).w(:,k)';          
        end;      
        J(sJ+1,n*(mm*ll+1)+1+2*n+nb_par+s1*(n+1)+s2*n+(1:s2))=...
          step_cnd(j).lambda_w';
        J(sJ+1,n*(mm*ll+1)+1+2*n+nb_par+s1*(n+1)+s2*(n+1)+(1:s1))=...
            (step_cnd(j).alpha)';      
        res(sJ+1)=0; 
      end;       
      % add extra conditions
      if method.extra_condition
        [resi,condi]=sys_cond(point);
        for j=1:length(condi)
          k=size(J,1);
          res(k+1)=resi(j);          
          for o=1:ll*mm+1
            J(k+1,(o-1)*n+1:o*n)=condi(j).profile(:,o)';
          end;        
          J(k+1,n*mm*ll+n+1)=condi(j).period;        
          for o=1:length(free_par)
            J(k+1,n*mm*ll+n+1+o)=condi(j).parameter(free_par(o));
          end; 
          nb_par=length(free_par);
          J(k+1,n*(mm*ll+1)+1+nb_par+(1:n))=condi(j).x1';
          J(k+1,n*(mm*ll+1)+1+nb_par+n+(1:n))=condi(j).x2';        
          s1=length(condi(j).v(1,:));
          for o=1:s1
            J(k+1,n*(mm*ll+1)+1+nb_par+2*n+(o-1)*n+(1:n))=condi(j).v(:,o)';
          end;        
          J(k+1,n*(mm*ll+1)+1+2*n+nb_par+s1*n+(1:s1))=condi(j).lambda_v';
          s2=size(condi(j).w,2);
          for o=1:s2
              J(k+1,n*(mm*ll+1)+1+nb_par+2*n+s1*(n+1)+(o-1)*n+(1:n))=...
                condi(j).w(:,o)';          
          end;        
          J(k+1,n*(mm*ll+1)+1+2*n+nb_par+s1*(n+1)+s2*n+(1:s2))=...
            condi(j).lambda_w';
          J(k+1,n*(mm*ll+1)+1+2*n+nb_par+s1*(n+1)+s2*(n+1)+(1:s1))=...
            condi(j).alpha';    
        end; 
      end;    
      % solve linear system
      if size(J,1)~=size(J,2)
        disp('P_CORREC warning: use of nonsquare Jacobian.');  
      end;
      dx=J\res;
      % apply non-parameter corrections
      for k=1:ll*mm+1
        point.profile(:,k)=point.profile(:,k)-real(dx((k-1)*n+1:k*n));
      end;
      point.period=point.period-real(dx(n*mm*ll+1+n));
      point.x1=point.x1-real(dx(n*mm*ll+1+n+(1:n)+length(free_par)));
      point.x2=point.x2-real(dx(n*mm*ll+1+2*n+(1:n)+length(free_par)));
      if ~isempty(point.v)
        s1=length(point.v(1,:));
      else 
        s1=0;
      end;
      for k=1:s1
        point.v(:,k)=point.v(:,k) ...
                   -dx(n*mm*ll+n+1+length(free_par)+2*n+(k-1)*n+(1:n));
        point.lambda_v(k)=point.lambda_v(k) ...
                   -dx(n*mm*ll+n+1+length(free_par)+2*n+s1*n+k);
      end;
      if ~isempty(point.w)
        s2=length(point.w(1,:));
      else 
        s2=0;
      end;
      for k=1:s2
        point.w(:,k)=point.w(:,k)-...
            dx(n*mm*ll+n+1+length(free_par)+2*n+(k-1)*n+(1:n)+s1*(n+1));
        % complex conjugate here because the Jacobian is formulated in terms of the complex 
        % conjugate of the residual... (NOT for the w's)
        point.lambda_w(k)=point.lambda_w(k)-conj(dx(n*mm*ll+n+1+length(free_par)+...
            2*n+s1*(n+1)+s2*n+k));  
      end;
      point.alpha=point.alpha-...
      dx(n*mm*ll+n+1+length(free_par)+2*n+length(point.v(1,:))*(n+1)+s2*(n+1)+(1:s1));
  end;

  % apply parameter corrections:
  for j=1:length(free_par)
    par(free_par(j))=par(free_par(j))-real(dx(p_start+j));
  end;

  % fill in parameters:
  point.parameter=par;  

  % compute residual:
  norm_res=norm(res);
  if print_r
    norm_residual=[i norm_res]
  end;

  % check for convergence
  if i==1
    r1=norm_res;
    r0=r1/2;
  else
    r0=r1;
    r1=norm_res;
  end;
  if r1<=conv_r | (i-1>nmon_iter & r1>=r0)
    break;
  end;

end;

% compute final residual

n_res=norm(res);

success=(n_res<=method.minimal_accuracy);

% recorrect with adapted mesh if necessary

if point.kind=='psol' | point.kind=='hcli'
  if size(point.mesh,2)
    ma=method.adapt_mesh_after_correct;
    if p_nr==1 | ( p_nr>1 & mod(p_nr,ma)==0 ) 
      % do not adapt mesh when p_nr=0
      % adapt mesh when p_nr=1
      % adapt mesh when p_nr>1 & mod(p_nr,ma)=0
      if success % adapt & correct again
        method2=method;
        method2.adapt_mesh_before_correct=1;
        method2.adapt_mesh_after_correct=0;
        [point,success]=p_correc(point,free_par,step_cnd,method2,2,previous);
      end
    end;
  end;
end;

return;
