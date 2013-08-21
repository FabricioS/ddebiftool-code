function root_plt(l0,l1,n1)

% function root_plt(l0,l1,n1)
% INPUT:
%	l0 approximated roots 
%	l1 corrected roots
%	n1 Newton convergence of l1 roots

% (c) DDE-BIFTOOL v. 1.02, 21/09/2001

hold on;

lr=[];
lg=[];
lk=[];

for i=1:length(l0),
  if real(l0(i))>0,
    lr(length(lr)+1)=l0(i); 
  else
    if real(l0(i))<0,
      lg(length(lg)+1)=l0(i);
    else
      lk(length(lk)+1)=l0(i);
    end;
  end;
end;

if length(lr),
  plot(real(lr),imag(lr),'rx');
end;
if length(lg),
  plot(real(lg),imag(lg),'gx');
end;
if length(lk),
  plot(real(lk),imag(lk),'kx');
end;

lr=[];
lg=[];
lk=[];

if isempty(n1),
  n1=ones(size(l1));
end;

for i=1:length(l1),
  if n1(i)~=-1,
    if real(l1(i))>0,
      lr(length(lr)+1)=l1(i); 
    else
      if real(l1(i))<0,
        lg(length(lg)+1)=l1(i); 
      else
        lk(length(lk)+1)=l1(i); 
      end;
    end;
  end;
end;

if length(lr),
  plot(real(lr),imag(lr),'r*');
end;
if length(lg),
  plot(real(lg),imag(lg),'g*');
end;
if length(lk),
  plot(real(lk),imag(lk),'k*');
end;

a=axis;
if a(1)<0 & a(2)>0,
  plot([0 0],[a(3) a(4)],'-.'); 
end;

plot([a(1) a(2)],[0 0],'-.');

return;
