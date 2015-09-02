function biftype = num2bif(index)
%% Convert an index to a bifurcation type
% Return number of supported types on 'count'
%
% $Id$
%
%%
if nargin>0
    biftype=bif_num(index,'<-');
else
    biftype=bif_num();
end
end
