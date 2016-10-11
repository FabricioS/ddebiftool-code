function stability=p_stabil(funcs,p,method)
%% compute stability information for point
% function stability=p_stabil(funcs,point,method)
% INPUT:
%   funcs problem functions
%	point solution point
%	method method parameters 
% OUTPUT:
%	stability stability information

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

if strcmp(p.kind, 'psol') || strcmp(p.kind, 'pdbl') || strcmp(p.kind, 'pfld')
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
  if ~funcs.tp_del
    mu=mult_app(funcs,p.period,p.profile,mesh,p.degree,rho,max_number,col,p.parameter);
  else 
    d_ac=method.delay_accuracy;
    mu=mult_app(funcs,p.period,p.profile,mesh,p.degree,rho,max_number,col,p.parameter,d_ac);
  end;
  if ~isempty(mu)
    [~, index_vector]=sort(abs(mu));
    stability.mu=mu(index_vector(length(index_vector):-1:1));
  else
    stability.mu=[];
  end;
else
  % Other kinds of points are processed using the method from Breda et al
  % (collocation scheme to approximate rightmost eigenvalues)
  if isfield(method,'psd_nodes')
    N = method.psd_nodes;
  else
    warning('P_STABIL: stability method field "psd_nodes" is not set, using default value of 100');
    N = 100;
  end;
  if length(sys_tau()) < 2
    stability = p_stabil_ndde(funcs,p,N);
  else
    error('p_stabil: Can not handle stability of NDDE with multiple delays.');
  end
end;

return;
