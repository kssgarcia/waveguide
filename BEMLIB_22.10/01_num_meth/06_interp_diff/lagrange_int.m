close all
clear all

%----
% data:
%----

%x=[-1:0.1:1];
%f=1./(1+25*x.^2);

%---
%N=5;
%x=[0.0 0.3 0.5 0.6 0.9 1.2];
%f=[0.0 0.4 0.4 0.3 0.7 0.5];
%f=[0.0 0.0 0.0 0.0 1.0 0.0];
%---

N=2;
step=2.0/N;
for i=1:N+1
 x(i)=-1.0+step*(i-1.0);
 f(i)= 1.0/(1.0+25.0*x(i)^2);
end

Nexact=128;
step=2.0/Nexact;
for i=1:Nexact+1
 xexact(i)=-1.0+step*(i-1.0);
 yexact(i)= 1.0/(1.0+25.0*xexact(i)^2);
end

%----
% prepare a graph:
%----

nplot=128;
step=(x(N+1)-x(1))/nplot;

for l=1:nplot+1
  xint(l)=x(1)+step*(l-1);
  yint(l)=lagrange(N,x,f,xint(l));
end

hold on
plot(x,f,'o','markersize',5,'color','red')
plot(xint,yint,'linewidth',1)
plot(xexact,yexact,':','linewidth',1)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis([-1.0 1.0 -1.5 1.5]);
box

%----
% done
%----


