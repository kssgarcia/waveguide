%close all;
clear all;
N=16;
N=2;
Ncl=N/2;
Ncl=N/4;
Ncl=N;
Ncl=2;

%---
% interpolation nodes
%---

for i=1:N+1
 t(i) = cos((i-0.5)*pi/(N+1));
% h(i)= 1.0/(1.0+25.0*t(i)^2);
 h(i)= exp(-2.0*t(i)^2);
end

%----
% Evaluate Ti(tj)
%----

for j=1:N+1
 TM(1,j)=1.0;
 TM(2,j)=t(j);
 for i=3:N+1
  TM(i,j)= 2.0*t(j)*TM(i-1,j)-TM(i-2,j);
 end
end

%----
% coefficients
%----

for i=1:N+1
 c(i)=0.0;
 for j=1:N+1
   c(i)=c(i)+h(j)*TM(i,j);
  end
 c(i)=2.0*c(i)/(N+1.0);
end
c(1)=0.5*c(1);

%----
% prepare a graph:
%----

nplot=128;
step=2.0/nplot;

for l=1:nplot+1
  tint(l) = -1.0+step*(l-1);
  yint(l) = cheb_clen(Ncl,c,tint(l));
end

hold on
%plot(t,h,'o','markersize',5,'color','red')
plot(tint,yint,':','linewidth',1)
%plot(tint,yint,'--','linewidth',1)
%plot(tint,yint,'linewidth',1)
xlabel('t','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis([-1.0 1.0 -1.5 1.5]);
axis([-1.0 1.0 -0.5 1.5]);
box
