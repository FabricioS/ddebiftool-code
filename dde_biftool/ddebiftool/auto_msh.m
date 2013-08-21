function [TMNEW]=auto_msh(NDIM,UPS,NOLD,NCOLD,TMOLD,DTMOLD,NNEW)

% function [tmnew]=auto_msh(ndim,ups,nold,ncold,tmold,dtmold,nnew)
% INPUT:
%	ndim system dimension
%	ups solution profile
%	nold old number of intervals
%	ncold number of collocation parameters
%	tmold old mesh
%	dtmold delta's of old mesh
%	nnew new number of mesh points
% OUTPUT:
%	tmnew new mesh
% COMMENT: 
%	this function is a matlab translation of the AUTO
%	fortran function NEWMSH restricted to the case of periodic 
%	boundary conditions

% (c) DDE-BIFTOOL v. 1.00, 15/03/2000

% put the values of the monotonely increasing function EQDF in EQF:

[EQF]=auto_eqd(NOLD,NDIM,NCOLD,DTMOLD,UPS);

% uniformly divide the range of EQDF:

NOLDP1=NOLD+1;
NNEWP1=NNEW+1;
DAL=EQF(NOLDP1)/NNEW;
for J=1:NNEWP1,
  UNEQ(J)=(J-1)*DAL;
end;

[IAL]=auto_ord(NOLDP1,EQF,NNEWP1,UNEQ);

% generate the new mesh in TMNEW:

for J1=1:NNEWP1,
  J=IAL(J1);
  X=(UNEQ(J1)-EQF(J))/(EQF(J+1)-EQF(J));
  TMNEW(J1)=(1.d0-X)*TMOLD(J)+X*TMOLD(J+1);
end;

TMNEW(NNEWP1)=1;

return;
