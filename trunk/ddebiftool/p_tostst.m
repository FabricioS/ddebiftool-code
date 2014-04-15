function [stst,stpcnd]=p_tostst(point,pert)

% function [stst_point,stpcnd]=p_tostst(point,pert)
% INPUT:
%	point 
%       % for branch switching from stst
%       pert size of the perturbation onto the emanating branch
% OUTPUT:
%	stst_point starting guess for stst point derived from point
% COMMENT:
%       if point kind is stst and the point is near a 
%         branching bifurcation, stst_point is a guess on the
%         branch and stpcnd should be used as an extra condition
%         in p_correc 
%       if point kind is fold or hopf, stst_point contains the
%         same point, but now as steady state only, no stpcnd
%         is returned
%       if point kind is hcli stst_point(1) is initial
%         and stst_point(2) is final steady state of the
%         given connecting orbit, no stpcnd is returned

% (c) DDE-BIFTOOL v. 2.03, 01/08/2002

stst.kind='stst';
stst.parameter=point.parameter;

switch point.kind,
  case 'stst',
    n=length(point.x); % system dimension
    m=length(sys_tau); % number of delays
    xx=point.x;
    for i=1:m
      xx=[xx point.x];
    end;
    res=sys_rhs(xx,point.parameter);
    J=zeros(n,n);
    for i=0:m
      J(1:n,1:n)=J(1:n,1:n)+sys_deri(xx,point.parameter,i,[],[]);
    end;
    [e1,e2]=eig(J);
    e=diag(e2);
    [i1,i2]=min(abs(e));
    if abs(imag(e(i2)))>0,
      disp('P_TOSTST: Warning: smallest eigenvalue is complex, taking real part!');
    end;
    stst.x=point.x+pert*real(e1(:,i2));
    stpcnd=stst;
    stpcnd.x=real(e1(:,i2));
    stpcnd.parameter=0*point.parameter; 
  case {'fold','hopf'}, 
    stst.x=point.x;
  case 'hcli', 
    stst(2)=stst(1);
    stst(1).x=point.x1;
    stst(2).x=point.x2; 
  case 'psol',
    error('P_TOSTST: conversion psol to stst not supported.');
  otherwise,
    err=point.kind,
    error('P_TOSTST: point kind not recognized.');
end;

return;
