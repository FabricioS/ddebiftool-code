function mult_plt(mu)

% function mult_plt(mu)
% INPUT:
%	mu approximated multipliers

% (c) DDE-BIFTOOL v. 1.00, 11/03/2000

theta=0:0.005:2*pi;

plot(cos(theta),sin(theta));
hold on;

rmu=[];
gmu=[];
kmu=[];

for i=1:length(mu)
  if abs(mu(i))>1
    rmu(length(rmu)+1)=mu(i); % plot(real(mu(i)),imag(mu(i)),'r*');
  else
    if abs(mu(i))<1
      gmu(length(gmu)+1)=mu(i); % plot(real(mu(i)),imag(mu(i)),'g*');
    else
      kmu(length(kmu)+1)=mu(i); % plot(real(mu(i)),imag(mu(i)),'k*');
    end;
  end;
end;

if length(rmu)
  plot(real(rmu),imag(rmu),'rx');
end;
if length(gmu)
  plot(real(gmu),imag(gmu),'gx');
end;
if length(kmu)
  plot(real(kmu),imag(kmu),'kx');
end;

return;