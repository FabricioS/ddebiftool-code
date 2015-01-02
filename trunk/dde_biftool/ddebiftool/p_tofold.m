function fold=p_tofold(funcs,point)
%% convert stst point to fold point
% function fold_point=p_tofold(funcs,point)
% INPUT:
%   funcs problem functions
%	point point near fold point
% OUTPUT:
%	fold_point starting guess for fold point derived from point

% (c) DDE-BIFTOOL v. 1.00, 08/04/2000
%
% $Id$
%
%%
fold.kind='fold';
fold.parameter=point.parameter;

switch point.kind
    case {'stst','hopf'}
        x=point.x;
    case 'fold'
        error('P_TOFOLD: point is already fold.');
    case 'psol'
        error('P_TOFOLD: periodic psol to fold not supported.');
    otherwise
        error('P_TOFOLD: point type %s to fold not supported.',point.kind);     
end
fold.x=x;
D=ch_matrix(funcs,x,point.parameter,0);
[E1,E2]=eig(D);
[i1,i2]=min(abs(diag(E2))); %#ok<ASGLU>
fold.v=real(E1(:,i2));
end
