function p_pplot(point,component,line_type)

% function p_pplot(point,component,colour)
% INPUT:
%       point point whose profile needs plotting
%	component optional component value
%	colour optional colour to plot with

% (c) DDE-BIFTOOL v. 1.02, 14/11/2000

if point.kind~='psol' & point.kind~='hcli'
  err=point.kind
  error('P_PPLOT: point does not contain a profile.');
end;

d=point.degree;
L=length(point.mesh);

if ~exist('line_type')
  line_type=' ';
end;

if ~exist('component');
  component=[];
end;

if ischar(component)
  line_type=component;
  component=[];
end;

if L==0
  L=size(point.profile,2);
  mesh=0:1/(L-1):1;
else
  mesh=point.mesh;
end;

if isempty(component) 
  plot(mesh,point.profile,line_type);
  hold on;
  plot(mesh(1:d:L),point.profile(:,1:d:L),strcat(line_type,'.'));
else
  plot(mesh,point.profile(component,:),line_type);
  hold on;
  plot(mesh(1:d:L),point.profile(component,1:d:L),strcat(line_type,'.'));
end;

return;



