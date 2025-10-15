%============================================
% plot the Fourier series of a saw-tooth wave
%============================================

clear all
close all

%---
% parameters
%---

L = 1.0;

Nf = 128;

alpha = 1.5;
alpha = 0.8;

%---
% number of plotting points
%---

Npl = 2*128;

%---
% plotting
%---

figure(1)
hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
axis([-2*L 2*L -1.4 1.4])
box on

figure(2)
hold on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
set(gca,'fontsize',15)
ylim = 0.4;
ylim = 0.05;
axis([-2*L 2*L -ylim ylim])
box on
plot([-2*L 2*L], [0,0],'k')
plot([-L -L], [-ylim,ylim],'k')
plot([ L  L], [-ylim,ylim],'k')

%---
% prepare
%---

Dx=2*L/Npl;

%---
% sum
%---

for ip=-1:2
for i=1:Npl+1
  x(i) = -L+ Dx*(i-1)+2*ip;
  y(i) = 0.0;
  y1(i) = 0.0;
  for p=1:Nf
    fc  = -(-1)^p*2/(p*pi);
    fc1 = -(-1)^p*2/(p*pi)*(L/(p*pi))^alpha;
    y(i) = y(i)+fc*sin(p*pi*x(i)/L);
    y1(i) = y1(i)+fc1*sin(p*pi*x(i)/L);
  end
end
figure(1)
plot(x,y,'k')
figure(2)
y1=y1/(2*pi);
plot(x,y1,'k')
end

xplt(1)=-L;
xplt(2)= L;
yplt(1)=-1;
yplt(2)= 1;
figure(1)
plot(xplt,yplt,'r','LineWidth',2)

