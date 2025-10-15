%----
clear all
close all
%----

%------
% prepare graphs of the error function
%------

N = 128;
a = 0.0;
b = 3.0;

Dx = (b-a)/N;

figure(1)
hold on

%---
for i=1:N+1
 x(i) = 0.001+ a +(i-1)*Dx;
 y(i) = erfun(x(i));
 z(i) = 1.0-y(i);
end
%---

%---
plot(x,y,'k');
plot(x,z,'r--');
%---

axis([0, 3, 0, 1.0])
xlabel('w','fontsize',15)
ylabel('erf(w)    erfc(w)','fontsize',15)
set(gca,'fontsize',15)
box
