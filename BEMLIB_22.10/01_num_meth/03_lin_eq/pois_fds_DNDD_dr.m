clear all
close all

%==========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%==========================================

%------------------------------------------------
% This program accompanies the book:
%             C. Pozrikidis
% Numerical Computation in Science and Engineering
%        Oxford University Press
%------------------------------------------------

%-----------------------------------------
% Driver for generating a linear system
% for the finite-difference solution
% of the Poisson equation lapl(f)+g = 0,
% as discussed in the text
%-----------------------------------------

%---
% input
%---

 ax = 0.0;
 bx = 1.0;

 ay = 0.0;
 by = 1.0;

 Nx = 16;
 Ny = 16;

%---
% prepare
%---

 Dx = (bx-ax)/Nx;
 Dy = (by-ay)/Ny;


%--------------------
% boundary conditions
%--------------------

  for j=1:Ny+1
    w(j) = 0.0D0;  % example
    q(j) =-5.0D0;  % example
  end

  for i=1:Nx+1
    z(i) = 0.5D0;  % example
    v(i) =-1.0D0;  % example
  end

%---
% source
%---

  for j=1:Ny+1
    for i=1:Nx+1
      g(i,j) = 0.0D0; % example
    end
  end

%---------------------------
% generate the linear system
%---------------------------

  [mats,mat,rhs] = pois_fds_DNDD  ...
             ...
   (ax,bx ...
   ,ay,by ...
   ,Nx,Ny ...
   ,g ...
   ,w,q,z,v ...
   );

%-----------------
% printing session
%-----------------

%    mat


%---
% solution
%---

  sol = rhs/mat';


%---
% distribute the solution
%---

 p = 0;     % counter

  for j=2:Ny
   for i=2:Nx+1
     p = p+1;
     f(i,j) = sol(p);
   end
  end

  for j=1:Ny+1
    f(   1,j) = w(j);
  end

  for i=1:Nx+1
    f(i,   1) = z(i);
    f(i,Ny+1) = v(i);
  end

%---
% grid nodes
%---

for j=1:Ny+1
   for i=1:Nx+1
     x(i,j) = ax + (i-1)*Dx;
     y(i,j) = ay + (j-1)*Dy;
   end
end

%---
% plot
%---

%mesh(x,y,f)
surf(x,y,f)
%axis([0 2 0 1 0 0.1])
axis equal
box on
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('f','fontsize',15)
set(gca,'fontsize',15)
view(-32,32)



