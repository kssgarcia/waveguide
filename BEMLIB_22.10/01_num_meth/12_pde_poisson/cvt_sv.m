close all;
clear all;

%=========================
% viscous flow in a rectangular cavity
% stream function/vorticity formulation
%=========================

%-----
% parameters
%-----

Ulid = 1.0;
Lx = 2.0;
Ly = 1.0;
Nx = 48;
Ny = 24;

visc = 1.005; % viscosity
visc = 0.005;
rho = 1.0;     % density

relax = 0.5;

Niteri = 5;       % number of inner iterations
Niterg = 1000;      % number of global iterations
tol = 0.0000001; % tolerance

vort_init=0.0;   % initial vorticity

%-------
% prepare
%--------

Dx = Lx/Nx;
Dy = Ly/Ny;
Dx2 = 2.0*Dx;
Dy2 = 2.0*Dy;
Dxs = Dx*Dx;
Dys = Dy*Dy;
beta = Dxs/Dys;
beta1= 2.0*(1.0+beta);

nu = visc/rho;  % kinematic viscosity

%-----------------------------------------
% initialize vorticity and stream function
%-----------------------------------------

for i=1:Nx+1
 for j=1:Ny+1
   x(i,j) = (i-1.0D0)*Dx;
   y(i,j) = (j-1.0D0)*Dy;
   psi(i,j) = 0.0;
   psi_new(i,j) = 0.0;
   vort(i,j) = -vort_init;
 end
end

%------------------
% global iterations
%------------------

for iter=1:Niterg

 save = vort;

%---------------------------------------
% Jacobi updating of the stream function
% at the interior nodes
%---------------------------------------

 for iteri=1:Niteri

  for j=2:Ny
   for i=2:Nx
    res = (psi(i+1,j)+psi(i-1,j)+ beta*psi(i,j+1)+beta*psi(i,j-1) ...
             +Dxs*vort(i,j))/beta1-psi(i,j);
    psi(i,j) = psi(i,j) + relax*res;
   end
  end

 end

%--------------------------------------------------
% Compute the vorticity at boundary grid points
% using the velocity boundary conditions
%
% Lower-order boundary conditions are commented out
%-------------------------------------

%---
% top and bottom walls
%---

 for i=2:Nx

%  vort(i,1  ) = 2.0*(psi(i,1)  -psi(i,2)) /Dys
%  vort(i,Ny+1) = 2.0*(psi(i,Ny+1)-psi(i,Ny))/Dys - 2.0*Ulid/Dys

   vort(i,1  ) = (7.0*psi(i,1)-8.0*psi(i,2) +psi(i,3))/(2.0*Dys);
   vort(i,Ny+1) = (7.0*psi(i,Ny+1)-8.0*psi(i,Ny)+psi(i,Ny-1))/(2.0*Dys) ...
                 - 3.0*Ulid/Dy;
  end

%---
% left and right walls
%---

 for j=2:Ny

%  vort(1  ,j) = 2.0*(psi(1,j)    - psi(2,j) )/Dxs
%  vort(Nx+1,j) = 2.0*(psi(Nx+1,j) - psi(Nx,j))/Dxs

  vort(1,j) = (7.0*psi(1,j)-8.0*psi(2,j)+psi(3,j))/(2.0*Dxs);
  vort(Nx+1,j) = (7.0*psi(Nx+1,j)-8.0*psi(Nx,j)+psi(Nx-1,j))/(2.0*Dxs);

 end

%--------------------------------
% compute the velocity at the grid points
%         by central differences
%--------------------------------

 for j=2:Ny
  for i=2:Nx
    ux(i,j) =   (psi(i,j+1)-psi(i,j-1))/Dy2;
    uy(i,j) = - (psi(i+1,j)-psi(i-1,j))/Dx2;
   end
 end

%--------------------------------------------
% iterate on Poisson's equation for the vorticity
%--------------------------------------------

 for iteri=1:Niteri

  for j=2:Ny
   for i=2:Nx
    source(i,j) =  ux(i,j)*(vort(i+1,j)-vort(i-1,j))/Dx2...
                 + uy(i,j)*(vort(i,j+1)-vort(i,j-1))/Dy2;
    source(i,j) = -source(i,j)/nu;
    res = (vort(i+1,j)+vort(i-1,j) + beta*vort(i,j+1)+beta*vort(i,j-1) ...
         +Dxs*source(i,j))/beta1-vort(i,j);
    vort(i,j) = vort(i,j) + relax*res;
   end
  end

 end        % of iteri
%--

%------------------
% monitor the error
%------------------

cormax=0.0;

for i=1:Nx+1
 for j=1:Ny+1
  res = abs(vort(i,j)-save(i,j));
  if(res>cormax)
    cormax = res;
  end
 end
end

cormax

if(cormax<tol)
 break
end

%=====
end % of iterg
%=====

%---
% % set up plotting vectors
%---

for i=1:Nx+1             
   xgr(i)=Dx*(i-1);
end
for j=1:Ny+1
   ygr(j)=Dy*(j-1);
end

%---
% vorticity plot
%---

surf(20*xgr,20*ygr,vort')
%contour(xgr,ygr,vort',32)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('\omega','fontsize',15)
set(gca,'fontsize',15)
axis([0 Lx 0 Ly -10 10])
axis equal

%---
% stream function contour plot
%---

figure
contour(xgr,ygr,psi',32)
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('\psi','fontsize',15)
set(gca,'fontsize',15)
axis([0 Lx 0 Ly])
axis equal
