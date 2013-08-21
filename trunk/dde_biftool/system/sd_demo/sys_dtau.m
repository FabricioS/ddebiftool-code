function dtau=sys_dtau(delay_nr,xx,par,nx,np)

% p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11

dtau=[];

% first order derivatives wrt state variables:
if length(nx)==1 & length(np)==0,
  if nx==0 % derivative wrt x(t)
    if delay_nr==3
      dtau(1:5)=0;
      dtau(2)=par(5)*par(10)*xx(2,2);
    elseif delay_nr==4
      dtau(1)=xx(2,3)/(1+xx(1,1)*xx(2,3))^2;
      dtau(2:5)=0;
    elseif delay_nr==5
      dtau(1:5)=0;
      dtau(4)=1;
    elseif delay_nr==6
      dtau(5)=1;
    else
      dtau(1:5)=0;
    end;
  elseif nx==1 % derivative wrt x(t-tau1)
    if delay_nr==3
      dtau(1:5)=0;
      dtau(2)=par(5)*par(10)*xx(2,1);
    else
      dtau(1:5)=0;
    end;
  elseif nx==2 % derivative wrt x(t-tau2)
    if delay_nr==4
      dtau(1:5)=0;
      dtau(2)=xx(1,1)/(1+xx(1,1)*xx(2,3))^2;
    else
      dtau(1:5)=0;
    end;
  else
    dtau(1:5)=0;
  end;
% first order derivatives wrt parameters:
elseif length(nx)==0 & length(np)==1,
  if delay_nr==1 & np==10
    dtau=1;
  elseif delay_nr==2 & np==11
    dtau=1;
  elseif delay_nr==3 & np==5
    dtau=par(10)*xx(2,1)*xx(2,2);
  elseif delay_nr==3 & np==10
    dtau=par(5)*xx(2,1)*xx(2,2);
  else 
    dtau=0;
  end;
% second order derivatives wrt state variables:
elseif length(nx)==2 & length(np)==0,
  dtau=zeros(5);
  if delay_nr==3
    if (nx(1)==0 & nx(2)==1) | (nx(1)==1 & nx(2)==0)  
      dtau(2,2)=par(5)*par(10);
    end;
  elseif delay_nr==4
    if nx(1)==0 & nx(2)==0
      dtau(1,1)=-2*xx(2,3)*xx(2,3)/(1+xx(1,1)*xx(2,2))^3;
    elseif nx(1)==0 & nx(2)==2   
      dtau(1,2)=(1-xx(1,1)*xx(2,3))/(1+xx(1,1)*xx(2,2))^3;
    elseif nx(1)==2 & nx(2)==0   
      dtau(2,1)=(1-xx(1,1)*xx(2,3))/(1+xx(1,1)*xx(2,2))^3;
    elseif nx(1)==2 & nx(2)==2   
      dtau(2,2)=-2*xx(1,1)*xx(1,1)/(1+xx(1,1)*xx(2,2))^3;
    end;
  end;
% mixed state parameter derivatives:
elseif length(nx)==1 & length(np)==1,
  dtau(1:5)=0;
  if delay_nr==3
    if nx==0 & np==5
      dtau(2)=par(10)*xx(2,2);
    elseif nx==0 & np==10 
      dtau(2)=par(5)*xx(2,2);
    elseif nx==1 & np==5
      dtau(2)=par(10)*xx(2,1);
    elseif nx==1 & np==10
      dtau(2)=par(5)*xx(2,1);
    end;
  end;
end;

if isempty(dtau)
  [delay_nr nx np]
  error('SYS_DTAU: requested derivative does not exist!');
end;

return;

