close all
clear all

%================
% simulation of random dart landing
%================

%---
% prepare
%---

figure(1)
hold on
xlabel('log(N)','fontsize',15)
ylabel('log(|\pi_{MC}-\pi|)','fontsize',15)
set(gca,'fontsize',15)
box

%---
% draw a line with slope -1/2
%---

xplot(1)= 0.0; yplot(1)= 0.0;
xplot(2)=10.0; yplot(2)=-5.0;

plot(xplot,yplot,'r');

%---
% Monte Carlo
%---

Nmax = 1024*2*2*2*2*2;

Ic=0;  % counter

jplot=32;

for N=1:Nmax

 x = 2.0*rand-1.0;
 y = 2.0*rand-1.0;

 if(x^2+y^2<1.0)
  Ic=Ic+1;
 end

 pi_mc=4.0*Ic/N;

 if(jplot==32)
  plot(log(N),log(abs(pi_mc-pi)),'.')
  jplot=0;
 end

 jplot=jplot+1;

end
