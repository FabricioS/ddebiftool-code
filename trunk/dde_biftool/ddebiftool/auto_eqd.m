function [EQF]=auto_eqd(NTST,NDIM,NCOL,DTM,UPS)

% function [eqf]=auto_eqd(ntst,ndim,ncol,dtm,ups)
% INPUT:
%	ntst number of intervals
%       ndim system dimension
%	ncol number of collocation points (per interval)
%       dtm delta's of mesh
%       ups solution profile
% OUTPUT:
%       eqf values of monotonically increasing function eqdf
% COMMENT: 
%       this function is a matlab translation of the AUTO
%       fortran function EQDF, except
%	- it is restricted to the case of periodic solutions
%       - it uses different ordering of the data in ups
%	- the integral is approximated using the trapezium rule

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

HMACH=1.0d-7;

[WH]=auto_cnt(NCOL);

SMALL=1;
for J=1:NTST,
  JP1=J*NCOL+1;
  SC=1/(DTM(J)^NCOL);
  for I=1:NDIM,
    HD(J,I)=WH(NCOL+1)*UPS(I,JP1);
    for K=1:NCOL,
      K1=K+(J-1)*NCOL;
      HD(J,I)=HD(J,I)+WH(K)*UPS(I,K1);
    end;
    HD(J,I)=SC*HD(J,I);
    if (abs(HD(J,I))>HMACH), SMALL=0; end;
  end;
end;

% Take care of "small derivative" case.

if (SMALL),
  for I=1:NTST+1,
    EQF(I)=I-1;
  end;
  return;
end;

for I=1:NDIM,
  HD(NTST+1,I)=HD(1,I);
end
DTM(NTST+1)=DTM(1);

% Compute approximation to (NCOL+1)-st derivative :

for J=1:NTST,
  JP1=J+1;
  DTAV=0.5*(DTM(J)+DTM(J+1));  
  SC=1/DTAV;
  for I=1:NDIM,
    HD(J,I)=SC*( HD(JP1,I)-HD(J,I) );
  end;
end;

% Define the equidistribution function :

PWR=1/(NCOL+1);
EQF(1)=0;
for J=1:NTST,
  EP=0;
  if J==1
    for I=1:NDIM,
      EP=EP+abs( HD(NTST,I) )^PWR;
    end;
  else
    for I=1:NDIM,
      EP=EP+abs( HD(J-1,I) )^PWR;
    end;
  end;
  E=0;
  for I=1:NDIM,
    E=E+abs( HD(J,I) )^PWR;
  end;
  EQF(J+1)=EQF(J)+DTM(J)*0.5*(E+EP);
end;

return;
