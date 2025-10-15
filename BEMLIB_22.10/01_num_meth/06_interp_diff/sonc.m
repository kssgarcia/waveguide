%---
% graph sonc
%---

close all
clear all

%---
N=9;
Nplt=64;
%---

x(1)=0.0;
x(N+1)=1.0;

L=x(N+1)-x(1);
h=L/N;

for i=1:N+1
 x(i) = x(1)+h*(i-1.0);
 ynull(i) = 0;
end

step=L/Nplt;

%---
% plot
%---

i=3;
i=2;

for k=1:Nplt+1
 xint=x(1)+step*(k-1.0)+0.001;
 w = xint-x(i);
 fsinc(k) = sin(pi*w/h)/(N*sin(pi*w/L));
 fsonc(k) = sin(pi*w/h)/(pi*w/h);
 lagr(k) = 1.0;
 for j=1:N+1
  if(j~=i)
   lagr(k)=lagr(k)*(xint-x(j))/(x(i)-x(j));
  end
 end
 xplot(k)=xint;
end

hold on
%plot(xplot,fsinc)
plot(xplot,fsonc)
plot(xplot,lagr,'--')
plot(x,ynull,'o','markersize',5,'color','red')
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
axis([0 1 -0.5 2.0])
set(gca,'fontsize',15)
box


