clear all
close all

%------------------------------
% prepare graphs of the
% complete elliptic integrals
%------------------------------

N=128;
a=0.0;
b=0.99;

Dx=(b-a)/N;

figure(1)
hold on

%---
 for i=1:N+1
  x(i) = 0.001+ a +(i-1)*Dx;
 [y(i),z(i)] = ell_int(x(i));
%[y1(i),z1(i)] = ellipke(x(i));
 end

 plot(x,y,'k');
 plot(x,z,'k--');
%plot(x,y,'r.');
%plot(x,z,'r.');
%---

%axis([0, 1, 0, 10])
%axis([0, 3, 0, 1.2])
xlabel('k','fontsize',15)
ylabel('F                  E','fontsize',15)
set(gca,'fontsize',15)
box
