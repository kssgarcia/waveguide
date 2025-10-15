clear all
close all

%==============
% draw a circle
%==============

Np = 64;

for i=1:Np+1
 angle = 2.0*pi*(i-1.0)/Np;
 x(i)  = cos(angle);
 y(i)  = sin(angle);
end

figure(1)
hold on
plot(x,y,'k');
plot(x,y,'r.');
xlabel('x','fontsize',13)
ylabel('y','fontsize',13)
set(gca,'fontsize',13)
axis square
box on

