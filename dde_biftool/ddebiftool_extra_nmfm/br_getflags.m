function spoints = br_getflags(branch)
%% Collect list of special point numbers in branch
%
% $Id$
%
%%
%#ok<*AGROW>
ll = length(branch.point);

spoints = [];
num = [];

for i = 1:ll
    f = branch.point(i).flag;
    if isempty(f) && ~strcmp(f,'')
       fprintf('Empty flag at point %d.\n', i);
    end
    ind = bif2num(f);
    if ind > 0
        if ind>length(num)
            num(ind)=0; 
        end
        num(ind) = num(ind) + 1;
        spoints(ind, num(ind)) = i;
    end
end

end

