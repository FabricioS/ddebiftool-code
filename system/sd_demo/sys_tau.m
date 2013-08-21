function tau=sys_tau(delay_nr,xx,par)

if delay_nr==1
  tau=par(10);
elseif delay_nr==2
  tau=par(11);
elseif delay_nr==3
  tau=2+par(5)*par(10)*xx(2,1)*xx(2,2);
elseif delay_nr==4
 tau=1-1/(1+xx(2,3)*xx(1,1));
elseif delay_nr==5
 tau=xx(4,1);
elseif delay_nr==6
 tau=xx(5,1);
end;

return;

