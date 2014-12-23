function index = bif2num(biftype)
%% Convert bifurcation type to its index
% Return number of supported types on 'count'
%
% $Id$
%%
if strcmp(biftype, 'stst')
   index = 0;
elseif strcmp(biftype,'hopf')
   index = 1;
elseif strcmp(biftype,'fold')
   index = 2;
elseif strcmp(biftype,'psol')
   index = 3;
elseif strcmp(biftype, 'hcli')
   index = 4;
elseif strcmp(biftype,'genh')
   index = 5;
elseif strcmp(biftype, 'hoho')
   index = 6;
elseif strcmp(biftype, 'zeho')
   index = 7;
elseif strcmp(biftype, 'count')
   index = 7;
elseif strcmp(biftype, '')
   index = -1;
else
   error('BIF2NUM: unknown bifurcation type "%s"!', biftype);
end

end

