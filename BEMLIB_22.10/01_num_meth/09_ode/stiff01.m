clear all
close all

%===
% integration of a stiff equation
% using RK1
%===

 t=0;
 x=0;
 tfinal=10;
 Dt = 0.02009;
 Nsteps = floor(tfinal/Dt)

 for i=1:Nsteps
  dxdt = fstiff(t,x);
  x = x + dxdt*Dt;
  t = t+Dt;
  tplot(i)=t;
  xplot(i)=x;
 end

hold on
%plot(tplot,xplot,'o')
plot(tplot,xplot,'-r')
xlabel('t','fontsize',15)
ylabel('x','fontsize',15)
set(gca,'fontsize',15)
axis([0 10 -1.2 1.2])

