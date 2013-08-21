function [xm,ym]=df_measr(stability,par_or_branch,kind)

% function [xm,ym]=df_measr(stability,par_or_branch,kind)
% INPUT:
%	stability nonzero for plotting stability information
%	par_or_branch free parameter list or branch
%	kind type of solution point
% OUTPUT:
%	xm default measure for x-axis
%	ym default measure for y-axis

% (c) DDE-BIFTOOL v. 2.00, 23/12/2000

if isfield(par_or_branch,'parameter')
  free_par=par_or_branch.parameter.free;
  kind=par_or_branch.point(1).kind;
else
  free_par=par_or_branch;
end;

if kind=='hcli' & stability
  error('DF_MEASR: heteroclinic solutions have no stability.');
end;

if length(free_par)==0
  if stability
    xm.field='stability';
    xm.row='all';
    xm.col=1;
    ym.field='stability';
    ym.row='all';
    ym.col=1;
    if kind=='psol'
      xm.subfield='mu';
      xm.func='real';
      ym.subfield='mu';
      ym.func='imag';
    else
      xm.subfield='l1';
      xm.func='real';
      ym.subfield='l1';
      ym.func='imag';
    end;
  elseif kind=='psol'
    xm.field='profile';
    xm.subfield='';
    xm.row=1;
    xm.col='all';
    xm.func='';
    ym.field='mesh';
    ym.subfield='';
    ym.row=1;
    ym.col='all';
    ym.func='';    
  elseif kind=='hcli'
    xm.field='profile';
    xm.subfield='';
    xm.row=1;
    xm.col='ampl';
    xm.func='';
    ym.field='period';
    ym.subfield='';
    ym.row=1;
    ym.col=1;
    ym.func='';  
  else
    xm.field='x';
    xm.subfield='';
    xm.row='all';
    xm.col=1;
    xm.func='';
    ym.field='x';
    ym.subfield='';
    ym.row='all';
    ym.col=1;
    ym.func='';  
  end;
else
  xm.field='parameter';
  xm.subfield='';
  xm.row=1;
  xm.col=free_par(1);
  xm.func='';
  if stability
    if kind=='psol'
      ym.field='stability';
      ym.subfield='mu';
      ym.row='all';
      ym.col=1;
      ym.func='abs';
    else
      ym.field='stability';
      ym.subfield='l1';
      ym.row='all';
      ym.col=1;
      ym.func='real';
    end;
  elseif length(free_par)>1
    ym.field='parameter';
    ym.subfield='';
    ym.row=1;
    ym.col=free_par(2);
    ym.func='';
  elseif kind=='psol'
    ym.field='profile';
    ym.subfield='';
    ym.row=1;
    ym.col='ampl';
    ym.func='';
  elseif kind=='hcli'
    ym.field='profile';
    ym.subfield='';
    ym.row=1;
    ym.col='ampl';
    ym.func='';
  else
    ym.field='x';
    ym.subfield='';
    ym.row=1;
    ym.col=1;
    ym.func='';
  end;
end;

return;
