function J=sys_dtau_SD_PObif(itau,x,p,nx,np,sys_dtau,dim,xtau_ind)
%% partial derivative of sys_tau for state-dependent delays
%
% $Id$
%
nvec=size(x,3);
xall=x(1:dim,:,:);
ntau=xtau_ind(1,end)-1;
if itau<=ntau
    J=sys_dtau(itau,xall(:,1:itau,:),p,nx,np);
else
    [ir,ic]=find(xtau_ind==itau+1,1,'first');
    %tau1=sys_tau(ir-1,xall(:,1:xtau_ind(ir,1)-1,:),p);
    J1=sys_dtau(ir-1,xall(:,1:xtau_ind(ir,1)-1,:),p,nx,np);
    %tau2=sys_tau(ic-1,xall(:,xtau_ind(ir,1:ic-1),:),p);
    J2=sys_dtau(ic-1,xall(:,xtau_ind(ir,1:ic-1),:),p,nx-(ir-1),np);
    %tau=tau1+tau2;
    J=J1+J2;
end
n_ext=size(x,1)-dim;
if length(nx)==1 && isempty(np)
    J=cat(2,J,zeros(1,n_ext,size(J,3)));
elseif length(np)==1 && isempty(nx)
    J=reshape(J,[1,1,nvec]);
elseif (length(nx)==2 && isempty(np)) || (length(nx)==1 && ~isempty(np))
    error('sys_dtau_SD_PObif:nx',...
        'higher derivatives of tau for PO bifurcations not implemented');
end
end
