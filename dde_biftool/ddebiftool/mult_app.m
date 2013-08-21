function mu=mult_app(period,profile,mesh,degree,rho,max_number,col,par,d_ac)

% function mu=mult_app(period,profile,mesh,degree,rho,max_number,col,par,d_ac)
% INPUT: 
%	period period of solution
%	profile periodic solution profile
%       mesh periodic solution mesh (if empty, mesh is assumed uniform)
%	degree degree of piecewise polynomials
%       rho keep multipliers with modulus >= rho
%	max_number keep at most max_number multipliers 
%	col collocation points
%       par current parameter values in R^p
%       d_ac (only for state-dependent delays) tau<d_ac is treated as 
%             tau<0 (stability is not computed)
% OUTPUT:
%       mu approximations of dominant multipliers

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

% if mesh is empty assume equidistant mesh:
if isempty(mesh)
  mesh=0:1/(size(profile,2)-1):1;
end;

m=degree;

tp_del=nargin('sys_tau');
if tp_del==0
  tau=par(sys_tau);
  d=length(tau);
else
  d=sys_ntau;
  tot_tau=[];
end;

if d>0
  if tp_del==0 
    tT=max(tau)/period;
  else
    % compute max(tau)
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

else

  ext_mesh=mesh;

end;

if tp_del==0
  M=mult_int(period,profile,mesh,m,ext_mesh,col,par);
else
  M=mult_int(period,profile,mesh,m,ext_mesh,col,par,d_ac);
end;

if isempty(M)
  mu=[];
  return;
end;

s=eig(M);

[dummy,I]=sort(abs(s));

mu=s(I(length(I):-1:max(1,length(I)-max_number+1)));

for i=length(mu):-1:1
  if abs(mu(i))>=rho
    break;
  end;
end;

mu=mu(1:i);

return;
