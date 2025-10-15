close all
clear all

%=======
% driver to draw stairs
%
% xs, ys, zs: left starting point
%
% steph: step height
% stepl: step length
% stepw: step width
%=======

xs = 0.0;
ys = 0.0;
zs = 0.0;

steph = 20.0;
stepl = 30.0;
stepw = 200.0;
nsteps = 10;

figure(1)

stairs (xs,ys,zs,steph,stepl,stepw,nsteps)

axis equal
view ([-24, 22])
xlabel('x')
ylabel('y')
zlabel('z')

box
