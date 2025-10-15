clear all;
close all;

%==================================
% conformal mapping of the exterior
% of a peridic array of cylinders
%==================================

alpha=0.6;
alpha=0.3;
alpha=0.8;
L=1.0;

Nxil=8;   % number of left xi divisions
Nxil=4;   % number of left xi divisions

Nxic=32;  % number of central xi divisions
Nxic=92;  % number of central xi divisions

Nxir=8;  % number of right xi divisions
Nxir=4;  % number of right xi divisions

Neta=32;  % number eta divisions

[xi,eta,x,y,h]= conf2_grid(alpha,L,Nxil,Nxic,Nxir,Neta);

Nxi=Nxil+Nxic+Nxir;  % total xi divisions

xi(Nxi+2)=xi(2)+L;   % wrap
eta(Nxi+2)=eta(2);


break

%===========================
% solve the Poisson equation
%===========================

%--------
% prepare
%--------

itermax=8000;
tol=0.00000001;
relax=1.9;

%===================================
% solve for Dirichlet top and bottom
%===================================

%---
% source
%---

for i=1:Nxi+1
 for j=1:Neta+1
  g(i,j)=1.0;
  g(i,j)=0.0;
 end
end

%---
% bottom and top boundary conditions
%---

for i=1:Nxi+1
 fbot(i)=0.0;
 ftop(i)=1.0;
 ftop(i)=0.0;
end

for i=Nxil+1:Nxil+Nxic+1
 fbot(i)=0.0;
 fbot(i)=1.0;
end

%---
% initialize
%---

for i=1:Nxi+2
 for j=1:Neta+1
  f(i,j) = 0.0;
 end
end

[f, iter, Iflag] = conf2_pois_gs_dpr ...
...
   (Nxi,Neta,xi,eta,x,y,h,g,itermax,tol,relax,fbot,ftop,f);

%---
% plot
%---

figure
hold on

for i=1:Nxi
 for j=1:Neta
 A=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)];
 B=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1)];
 C=[f(i,j),f(i+1,j),f(i+1,j+1),f(i,j+1)];
 patch(A,B,C,C);
 end
end
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('f','fontsize',15)
%axis equal
axis ([-0.5 0.5 0 1.5 0 1.0]);
box

%===============================
% solve for Dirichlet at bottom
% and homogeneous Neumann at top
%===============================

%---
% source
%---

for i=1:Nxi+1
 for j=1:Neta+1
  g(i,j)=1.0;
 end
end

%---
% bottom and top boundary conditions
%---

for i=1:Nxi+1
 fbot(i)=0.0;
 ftop(i)=0.0;
end

for i=Nxil+1:Nxil+Nxic+1
 fbot(i)=1.0;
 fbot(i)=0.0;
end

%---
% initialize
%---

for i=1:Nxi+2
 for j=1:Neta+1
  f(i,j) = 0.0;
 end
end

[f, iter, Iflag] = conf2_pois_gs_ndpr ...
...
   (Nxi,Neta,xi,eta,x,y,h,g,itermax,tol,relax,fbot,f);

%---
% plot
%---

figure
hold on

for i=1:Nxi
 for j=1:Neta
 A=[x(i,j),x(i+1,j),x(i+1,j+1),x(i,j+1)];
 B=[y(i,j),y(i+1,j),y(i+1,j+1),y(i,j+1)];
 C=[f(i,j),f(i+1,j),f(i+1,j+1),f(i,j+1)];
 patch(A,B,C,C);
 end
end
set(gca,'fontsize',15)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('f','fontsize',15)
%axis equal
%axis ([-0.5 0.5 0 1.5 0 0.2]);
box
