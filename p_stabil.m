function stability=p_stabil(p,method)

% function stability=p_stabil(point,method)
% INPUT:
%	point solution point
%	method method parameters 
% OUTPUT:
%	stability stability information

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

if p.kind=='psol' | p.kind=='pdbl' | p.kind=='pfld'
  if isempty(p.mesh)
    mesh=0:1/(size(p.profile,2)-1):1;
  else
    mesh=p.mesh;
  end;
  rho=method.minimal_modulus;
  max_number=method.max_number_of_eigenvalues;
  if isempty(method.collocation_parameters)
    col=poly_gau(p.degree);
  else
    col=method.collocation_parameters;
  end;
  if(nargin('sys_tau'))==0
    mu=mult_app(p.period,p.profile,mesh,p.degree,rho,max_number,col,p.parameter);
  else 
    d_ac=method.delay_accuracy;
    mu=mult_app(p.period,p.profile,mesh,p.degree,rho,max_number,col,p.parameter,d_ac);
  end;
  if length(mu)
    [y,index_vector]=sort(abs(mu));
    stability.mu=mu(index_vector(length(index_vector):-1:1));
  else
    stability.mu=[];
  end;
else
  if isfield(method,'psd_nodes')
    N = method.psd_nodes;
  else
    warning('P_STABIL: stability method field "psd_nodes" is not set, using default value of 100');
    N = 100;
  end;
  stability = p_stabil_ndde(p,N);
end;

return;
