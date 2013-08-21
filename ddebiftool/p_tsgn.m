function [delay_nr,tz]=p_tsgn(point)

% function [delay_nr,tz]=p_tsgn(point)
% INPUT:
%       point a point 
% OUTPUT:
%	delay_nr number of negative delay 
%       tz (for psol only) we want tau(tz)=0 and dtau/dt(tz)=0
 
% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

delay_nr=0;
tz=0;
d=sys_ntau; % number of delays 

if point.kind=='stst' | point.kind=='hopf' | point.kind=='fold' 
  x=point.x;
  xx=x;
  for j=1:d
    tau=sys_tau(j,xx,point.parameter);
    if tau<0
      delay_nr=j;
      break;
    end;
    xx=[xx x];
  end;
elseif point.kind=='psol'
  mm=point.degree;
  ll=(size(point.profile,2)-1)/mm;
  if isempty(point.mesh)
    mesh=0:1/(ll*mm):1;
  else
    mesh=point.mesh;
  end;
  % compute delays on period interval 
  for t_m=1:length(mesh)
    x=point.profile(:,t_m);
    xx=x;
    for j=1:d
      tau=sys_tau(j,xx,point.parameter);
      tau_all(j,t_m)=tau;
      tau=tau/point.period;
      t_tau=mesh(t_m)-tau;
      while t_tau<0,
        t_tau=t_tau+1;
      end;
      index_b=length(mesh)-mm;
      while (t_tau<mesh(index_b))
        index_b=index_b-mm;
      end;
      hhh_tau=mesh(index_b+mm)-mesh(index_b);
      t_tau_trans=(t_tau-mesh(index_b))/hhh_tau;
      Pb=poly_elg(mm,t_tau_trans);
      x_tau=point.profile(:,index_b:index_b+mm)*Pb';
      xx=[xx x_tau];
    end;
  end;
  % compute delays between representation points on period interval
  for i=1:length(mesh)-1
    mesh_new(i)=mesh(i)+0.5*(mesh(i+1)-mesh(i));
  end;
  for j=1:d
    tau_all_new(j,:)=psol_eva(tau_all(j,:),mesh,mesh_new,mm);
  end;
  % check sign of delay
  for j=1:d
    lrm=length(mesh);
    tau_all_ref(j,1:lrm)=tau_all(j,1:lrm);
    mesh_all_ref(j,1:lrm)=mesh;
    for i=1:length(mesh)-1
      if (tau_all_new(j,i)<tau_all(j,i)) & (tau_all_new(j,i)<tau_all(j,i+1))
        mesh_ref=mesh(i):(mesh(i+1)-mesh(i))/20:mesh(i+1);
        tau_all_ref(j,lrm+1:lrm+21)=spline(mesh,tau_all(j,:),mesh_ref);
        mesh_all_ref(j,lrm+1:lrm+21)=mesh_ref(1:21);
        lrm=lrm+21;
      end;
    end;
    [tau_n ind_tau]=min(tau_all_ref(j,:));
    if tau_n<0
      tz=mesh_all_ref(j,ind_tau);
      delay_nr=j;
      break;
    end;
  end;
elseif point.kind=='hcli'
  error('P_TSGN: this routine is not (yet) implemented for connecting orbits.'); 
else
  err=point.kind
  error('P_TSGN: point kind not recognized.');
end;

return;

