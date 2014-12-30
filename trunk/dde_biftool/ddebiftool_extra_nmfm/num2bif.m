function biftype = num2bif(index)
%% Convert an index to a bifurcation type
% Return number of supported types on 'count'
%
% $Id$
%
%%
pointtype={...
    'stst';... % steady state
    'hopf';... % Hopf point
    'fold';... % Fold of equilibria
    'psol';... % periodic orbit (incl local bifurcations)
    'hcli';... % connecting orbit
    'genh';... % generalized (degenerate) Hopf, Bautin bifurcation
    'hoho';... % Hopf-Hopf interaction (double Hopf bifurcation)
    'zeho'     % Fold-Hopf interaction, zero-Hopf, Gavrilov-Guckenheimer bifucation
    };
if nargin<1
    fprintf('known point types:\n');
    for i=1:length(pointtype)
        fprintf('type %d: %s\n',i-1,pointtype{i});
    end
    index='count';
end
if all(0<=index) && all(index<length(pointtype))
    if length(index)==1
        biftype=pointtype{index+1};
    else
        biftype=pointtype(index+1);
    end
elseif strcmp(index,'count')
      biftype = 7;
else
      error('NUM2BIF: unknown bifurcation index = %g', index);
end
end
