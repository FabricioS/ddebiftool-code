function spoints = br_getflags(branch)
% Give a list of special point numbers

ll = length(branch.point);

spoints = [];
num = zeros(1,10);

for i = 1:ll
    f = branch.point(i).flag;
    if isempty(f) && ~strcmp(f,'')
       fprintf('Empty flag at point %g.\n', i);
    end
    ind = bif2num(f);
    if ind > 0
        num(ind) = num(ind) + 1;
        spoints(ind, num(ind)) = i;
    end
end

end

