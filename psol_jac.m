function [J,res,tT,extmesh] = psol_jac(funcs, col, T, psol_prof,mesh,degree,par,free_par,phase,varargin)

% function [J,res]=psol_jac(c,T,psol_prof,t,deg,par,free_par,phase)
% INPUT:
%	c collocation parameters in [0,1]^deg
%	T period 
%	psol_prof profile in R^(n x deg*l+1)
%	t representation points in [0,1]^(deg*l+1)
%	deg degree piecewise polynomial
%	par current parameter values in R^p
%	free_par free parameters numbers in N^d 
%	phase use phase condition or not (s = 1 or 0)
%   wrapJ (optional key-value pair, default true) 
%        wrap time points periodically into [0,1]
% OUTPUT: 
%	J jacobian in R^(n*deg*l+n+s x n*deg*l+1+n+d)
%	res residual in R^(n*deg*l+n+s)
%   tT delays, scaled by period
%   extmesh mesh of time points, extended back to -max(tT(:))
%
% if wrapJ the returned jacobian is augmented with derivative wrt period and
% free_par and wrapped around periodically inside [0,1] 
% if ~wrapJ the returned jacobian puts its entries into the interval
%  [-min(delay,0),max(1-[0,delays])], no augmentation is done. This
%  Jacobian can be used for computation of Floquet multipliers and modes
%
% Heavily modified by David A.W. Barton (August 2005 to January 2006)
% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

default={'wrapJ',true};
options=dde_set_options(default,varargin,'pass_on');
sys_tau=funcs.sys_tau;
sys_ntau=funcs.sys_ntau;


% get the system dimension from the dimensions of profile
sysn = size(psol_prof,1);
% degree of interpolating polynomials
sysm = degree;
% number of mesh points (mesh points + representation points==length(mesh))
sysl = (length(mesh) - 1)/sysm;

% check that the mesh size is correct
if ((sysl ~= floor(sysl)) || (length(mesh) ~= size(psol_prof,2)))
    error('PSOL_JAC: Mesh size is wrong!');
end;

% check collocation points
if ((length(col) ~= 0) && (length(col) ~= sysm))
    error('PSOL_JAC: Wrong number of collocation points');
end;

% get number of delays to consider
tp_del=funcs.tp_del;
if (tp_del == 0)  
    % constant delays
    n_tau = sys_tau();    % location of delays in parameter list
    tau = par(n_tau);   % values of the delays
    tT = tau/T;         % normalised delay values
    nd = length(n_tau); % number of delays
else
    error('PSOL_JAC: State-dependent delay is not yet supported');
end;

% do some extra initialisation if we are working on an NDDE
% get the mass matrix M where M*x'(t) = f(x(t),x(t-tau),x'(t-tau))
if isfield(funcs, 'sys_ndde')
    massmatrix = sys_ndde();
    if ((size(massmatrix,1) ~= sysn) || (size(massmatrix,2) ~= sysn))
        if (size(massmatrix,1)*size(massmatrix,2)) == 1
            massmatrix = eye(sysn);
        else
            error('PSOL_JAC: Incorrectly sized mass matrix given by sys_ndde');
        end;
    end;
    ndde = 1;
else
    massmatrix = eye(sysn);
    ndde = 0;
end;

% create some collocation points (Gaussian) if we weren't given any
if (length(col) == 0)
    col = poly_gau(sysm);
    gauss_c = col;
    non_gauss = 0;
else
    gauss_c = poly_gau(sysm);
    non_gauss = 1;
end;

% phase condition initialisation
if phase
  gauss_abs=ones(1,sysm);
  g=poly_gau(sysm-1);
  for k=1:sysm
    for j=1:sysm-1
      gauss_abs(k)=gauss_abs(k)/(gauss_c(k)-g(j));
    end;
    for j=1:sysm
      if j~=k
        gauss_abs(k)=gauss_abs(k)/(gauss_c(k)-gauss_c(j));
      end;
    end;
  end;
  gauss_abs=gauss_abs/sum(gauss_abs);
end;

% create our Lagrange polynomials at the collocation points
P0 = zeros(sysm+1,sysm+1);
dP0 = zeros(sysm+1,sysm+1);
ddP0 = zeros(sysm+1,sysm+1);

for i = 1:sysm
    P0(i,:) = poly_elg(sysm,col(i));     % Lagrange polynomial
    dP0(i,:) = poly_del(sysm,col(i));    % and its 1st derivative
    ddP0(i,:) = poly_ddel(sysm,col(i));  % and its 2nd derivative
end;

% some optimisations
sysnm = sysn*sysm;
sysnml = sysnm*sysl;

% create the jacobian matrix and residual
J = zeros(sysnml + sysn + phase, sysnml + sysn + 1 + length(free_par));
% J = (collocation eqns + BCs + phase cond)x(interp pts + period + pars)
res = zeros(sysnml + sysn + phase,1);

% iterate through the mesh points
for l = 1:sysl

    % the starting index in the profile
    idx = sysm*(l - 1) + 1;
    % the starting index in the jacobian
    idxjac = (idx - 1)*sysn + 1;
    % find the start time
    ts = mesh(idx);
    % determine the width of the mesh interval
    h = mesh(idx + sysm) - ts;
    
    % iterate through the collocation points
    for m = 1:sysm
        
        % the starting index in the profile
        idxm = idx + m - 1;
        % the starting index in the jacobian
        idxjacm = idxjac + (m - 1)*sysn;
        % the range in res we want to modify (only a range if multidim)
        idxrange = [idxjacm:idxjacm + sysn - 1];
        % get the Lagrange polynomials on this interval
        P = P0(m,:);
        dP = dP0(m,:)/h;
        ddP = ddP0(m,:)/(h^2);
        % find the time of the current collocation point
        coltime = ts + col(m)*h; 
        % determine x(t) at the collocation point
        x = psol_prof(:,idx:(idx + sysm))*P';
        dx = psol_prof(:,idx:(idx + sysm))*dP';
        
        % phase_condition
        if phase && ~non_gauss
            fup = gauss_abs(m)*h*dx';
            i_l_i= (l - 1)*sysnm;
            for q=0:sysm
                qq=i_l_i+q*sysn;
                J(sysnml + sysn + 1,qq+1:qq+sysn) = J(sysnml + sysn + 1,qq+1:qq+sysn)+P(q+1)*fup;
            end;
        end;
    
        % determine x(t-tau) for all delays
        for numtau = 1:nd
        
            % where does t-tau look back to
            coltime_tau_i = mod(coltime - tT(numtau),1);
            % do a linear search for the corresponding mesh interval
            idxtau_i = length(mesh) - sysm;
            while (coltime_tau_i < mesh(idxtau_i))
                idxtau_i = idxtau_i - sysm;
            end;
            idxtau(numtau) = idxtau_i;
            % the width of the mesh interval
            h_tau = mesh(idxtau_i + sysm) - mesh(idxtau_i);
            % transform coltime_tau_i onto [0,1]
            coltime_trans = (coltime_tau_i - mesh(idxtau_i))/h_tau;
            % find the Lagrange polynomials
            Ptau(numtau,:) = poly_elg(sysm,coltime_trans);
            dPtau(numtau,:) = poly_del(sysm,coltime_trans)/h_tau;
            ddPtau(numtau,:) = poly_ddel(sysm,coltime_trans)/(h_tau^2);
            % compute x(t-tau)
            xtau(:,numtau) = psol_prof(:,idxtau_i:(idxtau_i+sysm))*Ptau(numtau,:)';
            dxtau(:,numtau) = psol_prof(:,idxtau_i:(idxtau_i+sysm))*dPtau(numtau,:)';
            ddxtau(:,numtau) = psol_prof(:,idxtau_i:(idxtau_i+sysm))*ddPtau(numtau,:)';
            
        % END: for numtau = 1:nd
        end;
        
        % NDDE: add derivatives for neutral equations (remember to rescale
        %   dxtau from [0:1] to [0:T])
        xx = [x xtau dxtau/T]; 
        % compute the RHS for the current collocation point
        f = sys_rhs(xx,par);
        % compute the residual for the collocation equations
        res(idxrange) = res(idxrange) + (massmatrix*dx - T*f);
        % add the free parameter derivatives to the jacobian
        for i = 1:length(free_par)
            df = sys_deri(xx,par,[],free_par(i),[]);
            J(idxrange,sysnml + sysn + 1 + i) = -T*df;
        end;
        % add (sum P'(c)*Delta u) to the jacobian
        for k = 0:sysm
            kk = idxjac + k*sysn;   
            J(idxrange,kk:(kk+sysn-1)) = J(idxrange,kk:(kk+sysn-1)) + massmatrix*dP(k+1);
        end;
        % add -f*Delta T to the jacobian
        J(idxrange,sysnml + sysn + 1) = J(idxrange,sysnml + sysn + 1) - f;
        % compute T*A0
        TA0 = T*sys_deri(xx,par,0,[],[]);
        % add (T*A_0*sum P(c)*Delta u) to the jacobian
        for k = 0:sysm
            kk = idxjac + k*sysn;
            J(idxrange,kk:(kk+sysn-1)) = J(idxrange,kk:(kk+sysn-1)) - TA0*P(k+1);
        end;
        
        % iterate through the delays and add to the jacobian
        for numdelay = 1:nd

            % compute A1
            A1 = sys_deri(xx,par,numdelay,[],[]);
            if ndde
                % NDDE: compute A2
                A2 = sys_deri(xx,par,numdelay + nd,[],[]);
            else
                A2 = 0;
            end;
            % determine where in the jacobian we need to modify
            idxtau_i = (idxtau(numdelay) - 1)*sysn;
            % add (-T*A1*sum P(\tilde{c})*Delta u) to the jacobian
            for k = 0:sysm
                kk = idxtau_i + k*sysn + 1;   
                J(idxrange,kk:(kk+sysn-1)) = J(idxrange,kk:(kk+sysn-1)) - A1*T*Ptau(numdelay,k+1) - A2*dPtau(numdelay,k+1);
            end;
            % add (-A1*tau*sum P'(\tilde{c})*u/T) to the jacobian
            J(idxrange,sysnml + sysn + 1) = J(idxrange,sysnml + sysn + 1) - A1*dxtau(:,numdelay)*tT(numdelay) - A2*(-dxtau(:,numdelay) + tT(numdelay)*ddxtau(:,numdelay))/T;
            
            % add extra derivatives in case the free parameter is the delay
            for i = 1:length(free_par)
                if (free_par(i) == n_tau(numdelay))
                    J(idxrange,sysnml + sysn + 1 + i) = J(idxrange,sysnml + sysn + 1 + i) + A1*dxtau(:,numdelay) + A2*ddxtau(:,numdelay)/T;
                end;
            end;
            
        % END: for numdelay = 1:nd
        end;
        
    % END: for m = 1:sysm
    end;

% END: for l = 1:sysl
end;

% periodicity conditions
J((sysnml + 1):(sysnml + sysn),1:sysn) = eye(sysn);
J((sysnml + 1):(sysnml + sysn),(sysnml + 1):(sysnml + sysn)) = -eye(sysn);
res((sysnml + 1):(sysnml + sysn)) = psol_prof(:,1) - psol_prof(:,end);

% phase conditions
if phase & non_gauss,
  for l_i=1:l
    index_a=(l_i-1)*sysm+1;
    for k=1:sysm
      fac=gauss_abs(k)*(mesh((l_i-1)*sysm+1)-mesh(l_i*sysm+1));
      dPa=poly_dla(mesh(index_a:index_a+sysm),gauss_c(k));
      Pa=poly_lgr(mesh(index_a:index_a+sysm),gauss_c(k));
      u_prime=psol_prof(:,index_a:index_a+m)*dPa';
      for q=1:sysm+1
        J(sysnml + sysn + 1,(l_i-1)*sysnm+1+(q-1)*sysn:(l_i-1)*sysnm+q*sysn)= ...
	  J(sysnml + sysn + 1,(l_i-1)*sysnm+1+(q-1)*sysn:(l_i-1)*sysnm+q*sysn) + fac*Pa(q)*u_prime';
      end;
    end;
  end;
end;
if phase
  res(sysnml + sysn + 1,1)=0;
end;


%FIXME
extmesh = mesh;
return;
