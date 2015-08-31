function nmfm_printbif_type(str,occurence)
%% shortcut printing out number of occurences of type str
%
% $Id$
%
if sum(occurence)>0
    fprintf('%d %s ',sum(occurence),str);
end
end
