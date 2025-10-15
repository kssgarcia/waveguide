clear all
close all

%=========================
% buffon needle simulation
%=========================

%--------
% settings
%--------

a=1.0;
b=0.9;
boa=b/a;

hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
axis square
set(gca,'fontsize',15)
box
axis([-0.5 1.5 -0.5 1.5])

%--------
% draw two horizontal lines
%--------

xplot(1) =-0.5*a;yplot(1) = 0.0;
xplot(2) = 1.5*a;yplot(2) = 0.0;
plot(xplot,yplot,'r','linewidth',2)

xplot(1) =-0.5*a;yplot(1) = a;
xplot(2) = 1.5*a;yplot(2) = a;
plot(xplot,yplot,'r','linewidth',2)

%--------
% Monte-Carlo
%--------

N=1024;
N=128;
N=32;
N=64;

M=0;

for i=1:N
 r1 = rand;
 r2 = rand;

 theta = pi*r2;
 cs = cos(theta);
 sn = sin(theta);

 if(r1<boa*sn)
  M = M+1;
 end

 %--------
 % Draw the needle:
 %--------

 xpos=    a*rand;  % random x position of the center
 ypos=0.5*a*r2;    % random y position of the center

 xplot(1) = xpos-0.5*b*cs;
 yplot(1) = ypos-0.5*b*sn;
 xplot(2) = xpos+0.5*b*cs;
 yplot(2) = ypos+0.5*b*sn;
 plot(xplot(1),yplot(1),'o');
 plot(xplot,yplot);
% yplot(1) = a-yplot(1);
% yplot(2) = a-yplot(2);
% plot(xplot(1),yplot(1),'o');
% plot(xplot,yplot);

end

P = M/N;
pi_mc = 2*boa/P
