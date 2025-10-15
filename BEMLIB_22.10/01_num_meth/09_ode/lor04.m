clear all
close all

%=================================
% integration of the Lorenz system
% of ODEs using RK4
%
% t: time
% Dt: time step
% Nsteps: Number of steps
%=================================

%---
% set parameters
%---

Dt = 0.01; 
Nsteps = 2^11;

%---
% initial condition
%---

x = [0.2 0.2 0.1]';

%---
% prepare
%---

figure(1)
hold on;
set(gca,'fontsize',15)
view(06,20);
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('z','fontsize',15)

%---
% initialize
%---

t = 0.0;

plot3 (x(1),x(2),x(3),'ko');
box on

%---
% run
%---

for k=1:Nsteps
 fk    = lor_fnc(t, x);
 xtmp1 = x + 0.5*Dt*fk;
 ftmp1 = lor_fnc(t+0.5*Dt, xtmp1);
 xtmp2 = x + 0.5*Dt*ftmp1;
 ftmp2 = lor_fnc(t+0.5*Dt, xtmp2);
 xtmp3 = x + Dt*ftmp2;
 ftmp3 = lor_fnc(t+Dt, xtmp3);
 ffina = (fk + 2.0*ftmp1 + 2.0*ftmp2 + ftmp3)/6.0;
 xnew  = x + Dt*ffina;
 plot3([x(1) xnew(1)] ,[x(2) xnew(2)] ,[x(3) xnew(3)],'k');
 x = xnew;
 t = t+Dt;
end 

plot3(x(1),x(2),x(3),'ro');

