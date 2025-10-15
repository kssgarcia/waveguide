%====================================
% FDLIB
%
% compute the molified color function
% in a 2D domain
%====================================

close all
clear all

%---
% parameters
%---

ax = 0.0; bx = 1.0;
ay = 0.0; by = 1.0;
Nx = 2*32; Ny = 2*32;
NSG = 64;

shape=1;  % circle
shape=2;  % sine

%---
% prepare
%---

L = bx-ax;
Dx = L/Nx;
Dy = (by-ay)/Ny;

%---
% Cartesian grid
%---

for i=1:Nx+1
 xgrid(i) = ax+(i-1)*Dx;
end

for j=1:Ny+1
 ygrid(j) = ay+(j-1)*Dy;
end

%---
% color function
%---

[clr Divcl Divcint] = color (ax,bx,ay,by,Nx,Ny,NSG,shape);

%---
% plot the color function
%---

figure(1)
hold on
zplot = clr(1:Nx+1,1:Ny+1);
mesh(xgrid,ygrid,zplot')
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('color','fontsize',15)
%axis([ax bx ay by -0.6 0.6])

%---
% another plot
%---

figure(2)
%zplot=Divcl(1:Nx+1,1:Ny+1);
zplot=Divcint(1:Nx+1,1:Ny+1);
mesh(xgrid,ygrid,zplot')
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('Div(color)','fontsize',15)

%===========================
% solve the Poisson equation
%===========================

for i=1:Nx+1
 for j=2:Ny
  g(i,j) = 0.0;  % zero source term
 end
end

%---
% conductivity
%---

for i=1:Nx+2
 for j=1:Ny+1
  k(i,j) = 0.01+1.0*clr(i,j);
 end
end

%---
% parameters
%---

itermax = 2000;
tol = 0.0000001;
relax = 1.0;

Tbot = 0.0;
Ttop = 1.0;

[T,iter,Iflag] = pois_gs_dprc ...
...
   (Nx,Ny,Dx,Dy,g,itermax,tol,relax,Tbot,Ttop,k);

%---
% plot the conductivity
%---

figure(3)
zplot=T(1:Nx+1,1:Ny+1);
mesh(xgrid,ygrid,zplot')
set(gca,'fontsize',14)
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
zlabel('T','fontsize',14)
%axis([ax bx ay by -0.6 0.6])
