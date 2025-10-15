clear all
%close all
hold on

%=====
% Recursive multigrid solution for the 1D Poisson eqn
%=====

a=0.0; b=1.0;
qL=0.0; fR=0.0;
ndiv=5;
nu1=3;
nu2=3;
ncycle=2;
omega=1/3;

L=b-a;
N=2^ndiv;
h=L/N;

%---
% initialize the fine-grid solution
%---

for i=1:N+1
 x(i) = a+(i-1)*h;
 f(i)=0.0;
 f(i)=0.1*sin(2*pi*x(i)/L);
 f(i)=0.1*rand-0.05;
end

f(N+1)=fR;

plot(x, f,'ro-');

%---
% graph
%---

hold on
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('f','fontsize',15)
axis([0 1 -0.1 1])
box

%---
% right-hand side of Af=b
%---

for i=1:N
 g(i)=exp(-2*x(i));
 g(i)=sin(2*pi*x(i)/L);
 g(i)=1.0;
 b(i)=h*h*g(i);
end
b(1)=b(1)-2*h*qL;
b(N)=b(N)+fR;

%---
% V cycles
%---

for cycle=1:ncycle
 f=mg_vcycle(N,f,b,nu1,nu2,omega);
 plot(x, f,'sc:');
end

plot(x, f,'sr:');
