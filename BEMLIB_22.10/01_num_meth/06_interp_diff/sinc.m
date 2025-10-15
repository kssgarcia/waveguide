%---
% graph sonc
%---

%----
close all
clear all
%----

N=11;
Nplt=256;

%---
% data
%---

x(1)=0.0;
x(N+1)=1.0;

L=x(N+1)-x(1);
h=L/N;

for i=1:N+1
 x(i) = x(1)+h*(i-1.0);
 ynull(i) = 0;
end

step=3*L/Nplt;

i=2;

%---
% plot
%---

for k=1:Nplt+1
 xint=-L+x(1)+step*(k-1.0)+0.001;
 w = xint-x(i);
 fsinc(k) = sin(pi*w/h)/(N*sin(pi*w/L));
 lagr(k) = 1.0;
 for j=1:N+1
  if(j~=i)
   lagr(k)=lagr(k)*(xint-x(j))/(x(i)-x(j));
  end
 end
 xplot(k)=xint;
end

hold on
plot(xplot,fsinc)
%plot(xplot,lagr,'--')
plot(x,ynull,'o','markersize',5,'color','red')
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
axis([-L 2*L -0.5 1.2])
set(gca,'fontsize',15)
box


