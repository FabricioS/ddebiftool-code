function funcs=set_funcs(varargin)
%% fill in funcs structure with user-defined functions for use with DDE-Biftool functions
%
% possible named arguments: 'sys_rhs' (mandatory),'sys_ntau', 'sys_cond',
% 'sys_deri', 'sys_dtau', 'x_vectorized' (logical)
%
% Examples
% for DDE x'=-x(t-tau)+a*x^3-x^5 with parameters [tau,a]:
% funcs=set_funcs(...
%   'sys_rhs',@(x,p)-x(1,2)+p(2)*x(1,1)^3-x(1,1)^5,...
%   'sys_tau',@()1)
% uses defaults df_deriv and df_derit for partial derivatives and no extra conditions
%
% $Id$
%
%% Process options
defaults={...
    'sys_rhs',[],...             % user-defined r.h.s.
    'sys_ntau',@()0,...          % number of delays (SD-DDEs only)
    'sys_tau',[],...             % index of delays as parameters or dependence of delays
    'sys_cond',@dummy_cond,...   % extra conditions
    'sys_deri',[],...            % Jacobians and second-order derivatives
    'sys_dtau',[],...            % Jacobians and second-order derivatives for delay (SD-DDEs only)
    'sys_mfderi',{},...          % higher-order derivatives (up to 5) for normal form computations
    'x_vectorized',false};
funcs=dde_set_options(defaults,varargin);
if isempty(funcs.sys_rhs)
    file=exist('sys_rhs','file');
    if file==2
        funcs.sys_rhs=@sys_rhs;
    else
        error('sys_rhs undefined');
    end
end
if isempty(funcs.sys_tau)
    file=exist('sys_tau','file');
    if file==2
        funcs.sys_tau=@sys_tau;
    else
        error('sys_tau undefined');
    end
end
%% test for state-dependent delay
funcs.tp_del=true;
try
    dummytau=funcs.sys_tau(); %#ok<NASGU>
    funcs.tp_del=false;
catch %#ok<CTCH>
    funcs.tp_del=true;
end
if funcs.x_vectorized
    funcs.sys_rhs=@(x,p)reshape(funcs.sys_rhs(x,p),[size(x,1),1,size(x,3)]);
    if funcs.tp_del
        funcs.sys_tau=@(itau,x,p)reshape(funcs.sys_tau(itau,x,p),[1,1,size(x,3)]);
    end        
end
if isempty(funcs.sys_deri)
    funcs.sys_deri_provided=false;
    funcs.sys_deri=@(x,p,nx,np,v)df_deriv(funcs,x,p,nx,np,v);
else
    funcs.sys_deri_provided=true;
    if funcs.x_vectorized
        funcs.sys_deri=@(x,p,nx,np,v)wrap_deri(x,p,nx,np,v,funcs.sys_deri);
    end
end
if isempty(funcs.sys_dtau)
    funcs.sys_dtau_provided=false;
    funcs.sys_dtau=@(it,x,p,nx,np)df_derit(funcs,it,x,p,nx,np);
else
    funcs.sys_dtau_provided=true;
    if funcs.x_vectorized
        funcs.sys_dtau=@(itau,x,p,nx,np)wrap_dtau(itau,x,p,nx,np,funcs.sys_dtau);
    end
end
if iscell(funcs.sys_mfderi)
    funcs.sys_mfderi_provided=false;
    mfargs=funcs.sys_mfderi;
    funcs.sys_mfderi=@(x,p,varargin)df_mfderiv(funcs,x,p,varargin,mfargs{:});    
else
    funcs.sys_mfderi_provided=true;    
end
end
function [resi,condi]=dummy_cond(point) %#ok
resi=[];
condi=[];
end
%% reshape for vectorization
function J=wrap_deri(x,p,nx,np,v,sys_deri)
n=size(x,1);
nvec=size(x,3);
J=sys_deri(x,p,nx,np,v);
J=reshape(J,[n,numel(J)/(n*nvec),nvec]);
end
function J=wrap_dtau(nr,x,p,nx,np,sys_dtau)
n=size(x,1);
nvec=size(x,3);
J=sys_dtau(nr,x,p,nx,np);
if length(nx)==1 && isempty(np)
    J=reshape(J,[1,n,nvec]);
elseif length(nx)==2
    J=reshape(J,[n,n,nvec]);
elseif  length(nx)==1 && length(np)==1
    J=reshape(J,[1,n,nvec]);
end
end
