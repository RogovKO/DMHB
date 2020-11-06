function eigenplot(eigenvalues)
% plot approximated roots and its corrections

%eigenvalue.l0 is the approximated roots 
%eigenvalue.l1 is the corrected roots by newton method
hold on;
lr=[];
lg=[];
lk=[];

for i=1:length(eigenvalues.l0),
  if real(eigenvalues.l0(i))>0,
    lr(length(lr)+1)=eigenvalues.l0(i); 
  else
    if real(eigenvalues.l0(i))<0,
      lg(length(lg)+1)=eigenvalues.l0(i);
    else
      lk(length(lk)+1)=eigenvalues.l0(i);
    end;
  end;
end;

if length(lr),
  plot(real(lr),imag(lr),'ro');
end;
if length(lg),
  plot(real(lg),imag(lg),'bo');
end;
if length(lk),
  plot(real(lk),imag(lk),'ko');
end;

lr=[];
lg=[];
lk=[];

for i=1:length(eigenvalues.l1),
    if real(eigenvalues.l1(i))>0,
      lr(length(lr)+1)=eigenvalues.l1(i); 
    else
      if real(eigenvalues.l1(i))<0,
        lg(length(lg)+1)=eigenvalues.l1(i); 
      else
        lk(length(lk)+1)=eigenvalues.l1(i); 
      end;
    end;
end;


if length(lr),
  plot(real(lr),imag(lr),'r+');
end;
if length(lg),
  plot(real(lg),imag(lg),'b+');
end;
if length(lk),
  plot(real(lk),imag(lk),'k+');
end;


a=axis;
if a(1)<0 & a(2)>0,
  plot([0 0],[a(3) a(4)],'k-.'); 
end;

plot([a(1) a(2)],[0 0],'k-.');
title('approximations (o), corrected root (+)', 'fontsize', 14)
return;




















