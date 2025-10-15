clear all
close all

%------
% prepare graphs of miscellaneous functions
%------

N=2*128;
a=0.0;
b=10.0;

Dx=(b-a)/N;

hold on

%---
for repeat=1:1

for i=1:N+1
 x(i) = a +(i-1)*Dx;
 xs = x(i)^2;
 if(repeat==1)
  t = x(i);
  y(i) = 1/1.0001 * (sin(t) - 0.01*cos(t)+ 0.01*exp(-100*t));
 elseif(repeat==2)
  y(i) = 1-xs;
 elseif(repeat==3)
  y(i) = 1-xs+xs*xs;
 elseif(repeat==4)
  y(i) = 1-xs+xs*xs-xs*xs*xs;
 elseif(repeat==5)
  y(i) = 1-xs+xs*xs-xs*xs*xs+xs*xs*xs*xs;
 end
end

if(repeat==1)
 plot(x,y);
 elseif(repeat==2)
 plot(x,y,'--');
 elseif(repeat==3)
 plot(x,y,':');
 elseif(repeat==4)
 plot(x,y,'--');
 elseif(repeat==5)
 plot(x,y,':');
end

end
%---

%axis([0, 1, 0, 10])
%axis([0, 3, 0, 1.2])
xlabel('t','fontsize',15)
ylabel('x','fontsize',15)
axis([0 10 -1.2 1.2])
set(gca,'fontsize',15)
box
