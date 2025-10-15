%close all
clear all
hold on

%========================
% integration of ODEs
% using different methods
%========================

%---
% definitions
%---

x0=1.0;

tot = 900.0;  % time interval
Dt=0.05;

Np = floor(tot/Dt)

x(1) = x0-x0^2*Dt;
time(1) = Dt;

%-----
% euler
%-----

for k=1:Np
 x(k+1) = x(k) + Dt*(k*Dt-x(k)^2);
 time(k+1) = time(k)+Dt;
end

%-----
% plotting
%-----

%plot(x0(1),x0(2),'s')
%plot(x(1,:),x(2,:),'-o')
%xlabel('x_1','fontsize',15)
%ylabel('x_2','fontsize',15)
%plot([0 time(1)],[x0(1),x(1,1)],'-o');plot(time,x(1,:),'-o')
%plot([0 time(1)],[x0,x(1)],'-o');plot(time,x,'o')
plot(time,x,'.')

xlabel('t','fontsize',15)
ylabel('x','fontsize',15)
%axis([0 10 -2 2])
set(gca,'fontsize',15)
