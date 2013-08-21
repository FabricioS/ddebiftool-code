function branch=df_brnch(free_par,kind)

% function branch=df_brnch(free_par,kind)
% INPUT:
%       free_par free parameter list
%       kind type of solution point
% OUTPUT:
%	branch empty branch with default values

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001

branch.method=df_mthod(kind);

branch.parameter.free=free_par;
branch.parameter.min_bound=[];
branch.parameter.max_bound=[];
branch.parameter.max_step=[];

tp_del=nargin('sys_tau');
if tp_del==0
  tau=sys_tau;
  for j=1:length(tau)
    branch.parameter.min_bound(j,:)=[tau(j) 0];
  end;
end;

branch.point=[];

return;
