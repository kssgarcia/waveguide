close all
clear all

%==========================
% random targeting of a disk
%===========================

%---
% options
%---

mode = 1;                     % uniform in the disk
mode = 3; sig1=1.0; sig2=0.5; % Gaussian 2D
mode = 2; mu=1.0;             % Guassian radial

%---
% prepare
%---

figure(1)
hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
axis equal
set(gca,'fontsize',15)
box

%---
% draw a circle
%---

for i=1:33
 theta = (i-1)*2*pi/32;
 xplot(i) = cos(theta);
 yplot(i) = sin(theta);
end

plot(xplot,yplot,'r','linewidth',2)

%---
% random points
%---

N=1024*2;

Ic=0;  % counter

for i=1:N

 if(mode==1)
   radius = sqrt(rand);
   theta = 2*pi*rand;
   x = radius*cos(theta);
   y = radius*sin(theta);
 elseif(mode==2)
   radius = -1/mu * log(rand);
   radius = sqrt(radius);
   theta = 2*pi*rand;
   x = radius*cos(theta);
   y = radius*sin(theta);
 elseif(mode==3)
   rad1 = sig1*sqrt(-2*log(rand));
   rad2 = sig2*sqrt(-2*log(rand));
   theta = 2*pi*rand;
   x = rad1*cos(theta);
   y = rad2*sin(theta);
 end

 figure(1)
 plot(x,y,'.');
end

axis equal
axis([-2 2 -2 2])
