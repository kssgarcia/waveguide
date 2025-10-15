clear all
close all
hold on

%===================================================
% Explicit multigrid solution for the 2D Poisson eqn
% Lf+g=0 in the square [a, b] x [a, b]
%
% with the homogeneous Dirichlet BC all around
%
% fine grid size is N=2^ndiv
%===================================================

a=0.0;
b=1.0;
ndiv=6;
nu1=3;
nu2=3;
ncycle=3;  % number of cycles

L=b-a;
N=2^ndiv;
h=L/N;
Nx=N;
Ny=N;

%---
% initialize the fine-grid solution
%---

for i=1:Nx+1
 x(i)=a+(i-1)*h;
end
for j=1:Ny+1
 y(j)=a+(j-1)*h;
end

f=zeros(Nx+1,Ny+1);

for j=2:Ny
 for i=2:Nx
   f(i,j)=0.0;
   f(i,j)=0.1*rand-0.05;
 end
end

%---
% graph
%---

%mesh(x,y,f);

hold on
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('f','fontsize',15)
axis([0 1 -0.1 1])
box

%---
% right-hand side of Af=b
%---

for j=1:Ny+1
 for i=1:Nx+1
  g(i,j)=exp(-2*x(i));
  g(i,j)=sin(2*pi*x(i)/L);
  g(i,j)=1.0;
  b(i,j)=h*h*g(i,j);
 end
end

%--------
% prepare
%--------

esave=zeros(ndiv,Nx+1,Nx+1);   % save the solution (f) and the error (e)
rsave=zeros(ndiv,Nx+1,Ny+1);   % save the residual (r)

%=========
% V cycles
%=========

for cycle=1:ncycle

%---
% presmoothing
%---

%nu1=200;
f = mg_gs(nu1,Nx,Ny,f,b);
%mesh(x,y,f);

%--
% residual
%--

 r=zeros(Nx+1,Ny+1);

 for j=2:Ny
  for i=2:Nx
   r(i,j)= f(i+1,j)+f(i-1,j)-4.0*f(i,j)+f(i,j+1)+f(i,j-1)+b(i,j);
  end
 end

 esave(1,:,:)=f;
 rsave(1,:,:)=r;

%----
% down to coarse
%----

 Nsys=N;

 for level=2:ndiv

  [xhalf, yhalf, rhalf] = mg_restrict(Nsys,Nsys,x,y,r);
  Nsys=Nsys/2;
  clear x r e;
  x=xhalf; y=yhalf; r=rhalf; e=zeros(Nsys+1,Nsys+1);

%---
  if(Nsys>2)
%---
  e = mg_gs(nu1,Nsys,Nsys,e,r);
%---
  else
%---
  e(2,2)=r(2,2)/4.0;
%---
  end
%---

  for j=2:Nsys
   for i=2:Nsys
    r(i,j)= e(i+1,j)+e(i-1,j)-4.0*e(i,j)+e(i,j+1)+e(i,j-1)+r(i,j);
   end
  end

  for i=1:Nsys+1
   for j=1:Nsys+1
    esave(level,i,j)=e(i,j);
    rsave(level,i,j)=r(i,j);
   end
  end

 end

%----
% up to fine
%----

  for level=ndiv-1:-1:1
   [xdouble, ydouble, edouble] = mg_prolongate(Nsys,Nsys,x,y,e);
   x=xdouble; y=ydouble;
   Nsys=2*Nsys;
   for k=1:Nsys+1
    for l=1:Nsys+1
     e(k,l) = esave(level,k,l)+edouble(k,l);
     r(k,l) = rsave(level,k,l);
   end
   end
  if(nu2>0)
   e = mg_gs(nu2,Nsys,Nsys,e,r);
   end
 end

 f=e;

%---
end
%---

 mesh(x,y,f);
