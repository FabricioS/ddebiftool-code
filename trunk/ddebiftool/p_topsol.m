function [psol,stpcond]=p_topsol(point,ampl,col_degree,nr_int)

% function [psol_point,stpcond]=p_topsol(point,ampl,col_degree,nr_int)
% INPUT:
%	point (with stability information if not hcli)
%	ampl amplitude of periodic solution guess
%       % from Hopf:
%	col_degree piecewise polynomial degree for periodic solution
%	nr_int number of intervals for mesh
%       % for period doubling:
%	col_degree collocation parameters (empty for Gauss) 
%       % from connecting orbit:
%       col_degree, nr_int are optional
% OUTPUT:
%	psol_point starting guess for periodic solution derived from point
%	stpcond steplength condition for use in correction of guess

% (c) DDE-BIFTOOL v. 2.02, 30/06/2002

psol.kind='psol';
psol.parameter=point.parameter;
stpcond.kind='psol';
stpcond.parameter=0*point.parameter;

switch point.kind,
  case 'hopf',
    mesh=0:1/(col_degree*nr_int):1;
    x=point.x;
    psol.mesh=mesh;
    psol.degree=col_degree;
    stpcond.mesh=mesh;
    stpcond.degree=col_degree;
    v=point.v/norm(point.v);
    for i=1:size(x,1)
      psol.profile(i,:)=x(i)+ampl*(real(v(i))*sin(2*pi*mesh)+imag(v(i))*cos(2*pi*mesh));
      stpcond.profile(i,:)=real(v(i))*sin(2*pi*mesh)+imag(v(i))*cos(2*pi*mesh);
    end;
    if abs(point.omega)>0
      psol.period=abs(2*pi/point.omega);
    else
      disp('P_TOPSOL: zero frequency in Hopf point, period set to zero.');
    end;
    stpcond.period=0;
  case 'psol',
    % remove trivial multiplier
    [i1,i2]=min(abs(point.stability.mu-1));
    point.stability.mu(i2) = 0;
    % find near-bifurcation multiplier
    [i1,i2]=min(abs(abs(point.stability.mu)-1));
    if real(point.stability.mu(i2))<0 & abs(point.stability.mu(i2)+1)<0.3
      [eig_mesh,eig_v]=mult_dbl(point.period,point.profile, ...
          point.mesh,point.degree,col_degree,point.parameter);
      psol.mesh=eig_mesh;
      psol.degree=point.degree;
      psol.profile=point.profile;
      ll=size(psol.profile,2);
      psol.profile(:,ll+1:2*ll-1)=point.profile(:,2:ll);
      psol.profile=psol.profile+ampl*eig_v;
      psol.period=2*point.period;
    elseif real(point.stability.mu(i2))>0 & abs(point.stability.mu(i2)-1)<0.3
      [eig_mesh,eig_v]=mult_one(point.period,point.profile, ...
          point.mesh,point.degree,col_degree,point.parameter);
      psol.mesh=eig_mesh;
      psol.degree=point.degree;
      psol.profile=point.profile+ampl*eig_v;
      psol.period=point.period;
    else
      error('P_TOPSOL: periodic solution is not close enough to bifurcation.');
      return;
    end;
    stpcond.mesh=eig_mesh;
    stpcond.degree=point.degree;
    stpcond.profile=eig_v;
    stpcond.period=0;
  case 'hcli',
    psol.mesh=point.mesh;
    psol.degree=point.degree;
    psol.profile=point.profile;
    psol.period=point.period;
    stpcond.mesh=point.mesh;
    stpcond.degree=point.degree;
    stpcond.profile=zeros(size(point.profile));
    stpcond.period=1;
  case 'stst', 
    error('P_TOPSOL: stst to psol not supported, convert to hopf first.');
  case 'fold',
    error('P_TOPSOL: fold to psol not supported, convert to hopf first.');
  otherwise,
    err=point.kind
    error('P_TOPSOL: point kind not recognized.');
end;

return;

