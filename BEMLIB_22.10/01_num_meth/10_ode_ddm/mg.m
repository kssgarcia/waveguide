clear all
close all
hold on

%===================================================
% Explicit multigrid solution for the 1D Poisson eqn
% in the interval [a, b]
%
% with Neumann BC at the left end: g'=qL
% and Dirichlet BC at the right end: f=fR
%
% fine grid size is N=2^ndiv
%===================================================

a=0.0;
b=1.0;
qL=0.0;
fR=0.0;
ndiv=6;
nu1=3;
nu2=3;
ncycle=3;  % number of cycles
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
 f(i)=0.5*rand-0.05;    % random
end

f(N+1)=fR;

%---
% graph
%---

plot(x, f,'ro-');

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

%--------
% prepare
%--------

esave=zeros(ndiv,N+1);   % save the solution (f) and the error (e)
rsave=zeros(ndiv,N);     % save the residual (r)

%=========
% V cycles
%=========

for cycle=1:ncycle

%---
% presmoothing
%---

f = mg_rich(nu1,omega,N,f,b);

%--
% residual
%--

 r(1)= b(1)-2.0*f(1)+2.0*f(2);
 for i=2:N
   r(i)= b(i)+f(i-1)-2.0*f(i)+f(i+1);
 end

esave(1,:)=f;
rsave(1,:)=r;

%----
% down to coarse
%----

 Nsys=N;

 for level=2:ndiv

  [xhalf, rhalf] = mg_restrict(Nsys,x,r);
  Nsys=Nsys/2;
  clear x r e;
  x=xhalf; r=rhalf; e=zeros(Nsys+1,1);

%---
  if(Nsys>2)
     e = mg_rich(nu1,omega,Nsys,e,r);
  else
     e(1) = r(1)+r(2);
     e(2) = 0.5*(r(1)+2.0*r(2));
     e(3)=0.0;
  end

%---
% update the residual
%---

   r(1)= r(1)-2.0*e(1)+2.0*e(2);
   for i=2:Nsys
     r(i)= r(i)+e(i-1)-2.0*e(i)+e(i+1);
   end

   for k=1:Nsys
    esave(level,k)=e(k);
    rsave(level,k)=r(k);
   end
%  plot(x, e,'o:');

 end

%----
% up to fine
%----

 for level=ndiv-1:-1:1
  [xdouble, edouble] = mg_prolongate(Nsys,x,e);
  x=xdouble;
  Nsys=2*Nsys;
  for k=1:Nsys
   e(k) = esave(level,k)+edouble(k);
   r(k) = rsave(level,k);
  end
  e(Nsys+1)=0.0;
  if(nu2>0)
   e = mg_rich(nu2,omega,Nsys,e,r);
   end
 end

 f=e;

 plot(x, f,'ro:');

%---
end
%---

%plot(x, f,'ro:');
