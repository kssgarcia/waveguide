close all
clear all

%===========================
% simulation of dart landing
%===========================

%---
% prepare
%---

figure(1)
hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
axis square
set(gca,'fontsize',15)
box

%---
% draw a circle
%---

ncrc = 32;

for i=1:ncrc+1
 theta = (i-1)*2*pi/ncrc;
 xplot(i) = cos(theta);
 yplot(i) = sin(theta);
end

plot(xplot,yplot,'r','linewidth',2)

%---
% Monte Carlo
%---

N = 1024*2*2;

Ic = 0;  % counter of inside landings

for i=1:N
 x = 2.0*rand-1.0;
 y = 2.0*rand-1.0;
 plot(x,y,'.');
 if(x^2+y^2<1)
  Ic = Ic+1;
 end
end

%---
% estimated number pi
%---

pi_mc = 4.0*Ic/N
