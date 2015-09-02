function index = bif2num(biftype)
%% Convert bifurcation type to its index
% Return number of supported types on 'count'
%
% $Id$
%%
if nargin>0
    index=bif_num(biftype,'->');
else
    index=bif_num();
end
end
