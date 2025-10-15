close all
clear all

%----
% parameters
%----

figure(1)

a=1.0;   % particle radius

rho1=0.0;
rho2=1.0;
kappa =-0.5;
kappa = 0.0;
kappa = 0.5;
kappa = 0.9;

smax=6.0;
ndiv=128;
ndiv=64;

%==============================

for Iloop=1:4

if(Iloop==1)
 capl=0.7 % capilary length
elseif(Iloop==2)
 capl=1.0 % capilary length
elseif(Iloop==3)
 capl=2.0 % capilary length
elseif(Iloop==4)
 capl=5.0 % capilary length
end

alpha=pi/8;
alpha=0.98*pi;
alpha=0.00;    % to unlock
alpha=0.02*pi;
alpha=999;    % to unlock

alphain=alpha;

[s,x,q,beta,alpha,xc,xcl,scl,al_plot,xc_plot]  ...
...
= flsphere(a,capl,rho1,rho2,kappa,alphain,smax,ndiv);

if(Iloop==1)
 al_plot1=al_plot;
 xc_plot1=xc_plot;
elseif(Iloop==2)
 al_plot2=al_plot;
 xc_plot2=xc_plot;
elseif(Iloop==3)
 al_plot3=al_plot;
 xc_plot3=xc_plot;
elseif(Iloop==4)
 al_plot4=al_plot;
 xc_plot4=xc_plot;
end

%=========
% plot the hydrostatic shape
%=========

Iskip=1;
if(Iskip==0)

close all
hold on
plot(x,s,'o')

plot([-2*a 2*a],[0 0])
plot([0 0],[-6*a 6*a])

axis equal
axis([-2 2 -2 2])
axis([-6 6 -6 6])
set(gca,'fontsize',15)
xlabel('y/a','fontsize',15)
ylabel('x/a','fontsize',15)

end

end

%=========
% plot xc vs alpha
%=========

figure(1)
hold on
plot(al_plot1,xc_plot1)
plot(al_plot2,xc_plot2)
plot(al_plot3,xc_plot3)
plot(al_plot4,xc_plot4,':')
set(gca,'fontsize',15)
xlabel('\alpha/\pi','fontsize',15)
ylabel('x_c/a','fontsize',15)
