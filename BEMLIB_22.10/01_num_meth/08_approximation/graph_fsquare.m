clear all
close all

%==============================================
% plot the sine Fourier series of a square wave
%==============================================

a = -pi;
b = 3*pi;

Nf = 65

Npl = 2*128;

%---
% prepare
%---

Dx =(b-a)/Npl;

%---
% plotting
%---

figure(1)
hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis([-2 8 -1.4 1.4])
box on

for i=1:Npl+1
  x(i) = a+ Dx*(i-1);
  y(i) = 0.0;
  for j=1:2:Nf
    y(i) = y(i)+sin(j*x(i))/j;
  end
  y(i) = 4.0*y(i)/pi;
end

plot(x,y)

xplt(1)=a;
yplt(1)=0;
xplt(2)=0;
yplt(2)=0;
xplt(3)=0;
yplt(3)=1;
xplt(4)=pi;
yplt(4)=1;
xplt(5)=pi;
yplt(5)=-1;
xplt(6)=2*pi;
yplt(6)=-1;
xplt(7)=2*pi;
yplt(7)=0;
xplt(8)=3*pi;
yplt(8)=0;
plot(xplt,yplt,'LineWidth',2)

