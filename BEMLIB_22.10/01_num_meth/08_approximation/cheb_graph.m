clear all
close all

%------
% prepare graphs of functions
% cheb poly are implemented
%------

N=128;
a=0.0;
b=0.99;

Dx=(b-a)/N;

hold on

%---
for repeat=1:5

for i=1:N+1
 x(i) = 0.001+ a +(i-1)*Dx;
 xs = x(i)^2;
 if(repeat==1)
  y(i) = 1/(1+xs);
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
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
box
