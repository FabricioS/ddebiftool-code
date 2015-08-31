function index = bif2num(biftype)
%% Convert bifurcation type to its index
% Return number of supported types on 'count'
%
% $Id$
%%
biflist={'stst','hopf','fold','psol','hcli','genh','hoho','zeho','BT','CP'};
index=find(strcmp(biftype,biflist))-1;
if isempty(index)
    switch biftype
        case 'count'
            index=length(biflist)-1;
        case ''
            index=-1;
        otherwise
            error('BIF2NUM: unknown bifurcation type "%s"!', biftype);
    end
end
end

