close all
clear all

%---
% parameters
%---

ax=0; bx=1.4;
ay=0; by=0.8;
Nx=32; Ny=32;

itermax=128*128;
itermax=12*12*128;

tol=0.00000001;
relax=1.5;
relax=1.0;
relax=0.5;


Dx=(bx-ax)/Nx;
Dy=(by-ay)/Ny;

%---
% source
%---

for i=1:Nx+1
 for j=1:Ny+1
   g(i,j)=0.0;
   g(i,j)=1.0;
   g(i,j)=sin(2*pi*(i-1)/Nx)*sin(2*pi*(j-1)/Ny);
 end
end

%---
% boundary conditions
%---

qbot= 1.0;
qtop=-1.0;
qleft=1.0;
qright=-1.0;

qleft=0.0;
qright=0.0;
qbot=0.0;
qtop=0.0;

%---
% initialize
%---

for i=1:Nx+1
 for j=1:Ny+1
  f(i,j) = 0.0;
 end
end

%---
% solve
%---

Ishift=1;

[f, iter, Iflag] = pois_gs_nnnn1 ...
 (Nx,Ny,Dx,Dy,g,itermax,tol,relax ...
 ,qleft,qright,qbot,qtop,f,Ishift);

%---
% plotting
%---

for i=1:Nx+1
 x(i)=ax+(i-1)*Dx;
end

for j=1:Ny+1
 y(j)=ay+(j-1)*Dy;
end

mesh(x,y,f') 
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('f','fontsize',15)
axis([ax bx ay by])
