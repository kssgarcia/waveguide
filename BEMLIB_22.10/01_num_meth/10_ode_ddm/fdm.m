clear all
close all

%---
% parameters
%---

N=32;
h=0.5*pi/N;
eps=10.0;

%---
% linear system
%---

for i=1:N
  t(i) = (i-1.0)*h;
  b(i) = 1.0;
  c(i) = 1.0;
  a(i) = -2.0+h*h*(4.0+eps*cos(t(i)));
  s(i) = 0.0;
end

b(1)=2.0;
s(1) = -4.0*h;
s(N) = 1.0;

x = thomas (N,a,b,c,s);

t(N+1)=0.5*pi;
x(N+1)=-1;

%-----
% plot
%-----

plot(t, x,'-o');
xlabel('x','fontsize',15);
ylabel('f','fontsize',15);
set(gca,'fontsize',15)
axis([0 0.5*pi -4 4])
