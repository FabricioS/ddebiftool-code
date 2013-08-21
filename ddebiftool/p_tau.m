function [tau_eva]=p_tau(point,d_nr,t)

% function [tau_eva]=p_tau(point,d_nr,t)
% INPUT:
%       point a point
%       d_nr number(s) of delay(s) to evaluate
%       t (optional) point(s) where to evaluate
% OUTPUT:
%       tau_eva value(s) of evaluated delay(s)

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

d_nr_eva=length(d_nr); % number of delays to evaluate
max_d_nr=max(d_nr); % maximum number of delay to evaluate
tau_eva=[];

if point.kind=='stst' | point.kind=='hopf' | point.kind=='fold' 
  x=point.x;
  xx=x;
  for jj=1:max_d_nr
    tau=sys_tau(jj,xx,point.parameter);
    for ii=1:d_nr_eva
      if jj==d_nr(ii)
        tau_eva(ii)=tau;
      end;
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
  if exist('t')
    mesh_eva=t;
    l_me=length(t);
  else
    mesh_eva=mesh;
    l_me=length(mesh);
  end;
  % compute delays 
  for t_m=1:l_me
    % compute xx at mesh_eva(t_m)
    xx=psol_eva(point.profile,mesh,mesh_eva(t_m),mm);
    % compute delays at mesh_eva(t_m)
    for jj=1:max_d_nr  
      tau=sys_tau(jj,xx,point.parameter);
      for ii=1:d_nr_eva
        if jj==d_nr(ii)
          tau_eva(ii,t_m)=tau;
        end;
      end;
      tau=tau/point.period;
      t_tau=mesh_eva(t_m)-tau;
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
elseif point.kind=='hcli'
  error('P_TAU: this routine is not (yet) implemented for connecting orbits.'); 
else
  err=point.kind
  error('P_TAU: point kind not recognized.');
end;

return;


























