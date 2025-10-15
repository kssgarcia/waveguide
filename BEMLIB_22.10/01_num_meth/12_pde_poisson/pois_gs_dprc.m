function [f, iter, Iflag] = pois_gs_dprc ...
...
   (Nx,Ny,Dx,Dy,g,itermax,tol,relax,fbot,ftop,k)

%------------------------------------------
% Solution of a modified Poisson equation
% in a rectangular domain
% with a uniform Dirichlet boundary condition
% at the top and bottom boundaries
% and periodic boundary condition
% at the left and right boundaries
%
% The solution is found by point Gauss-Seidel
% iterations
%------------------------------------------

%--------
% prepare
%--------

Dxs = Dx*Dx;
Dys = Dy*Dy;
beta = Dxs/Dys;
beta1 = 2.0*(1.0+beta);

Iflag = 0;  % convergence flag; 0 is good

%---
% initialize
%---

for i=1:Nx+2
 for j=1:Ny+1
  f(i,j) = 0.0;
 end
end

%---
% implement the Dirichlet
% boundary condition
%---

for i=1:Nx+2
 f(i, 1)   = fbot;  % bottom
 f(i,Ny+1) = ftop;  % top
end

%------------------------
% Gauss-Siedel iterations
%------------------------

for iter=1:itermax

%------------------------
% update nodes row-by-row
%------------------------

%---
% save
%---

 fsave = f;

%---
% iterate
%---

 cormax=0.0;

 for j=2:Ny
  for i=2:Nx+1
   kl=0.5*(k(i-1,j)+k(i,j));
   kr=0.5*(k(i+1,j)+k(i,j));
   kb=0.5*(k(i,j-1)+k(i,j));
   kt=0.5*(k(i,j+1)+k(i,j));
   res = ( kr*f(i+1,j) ...
          +kl*f(i-1,j) ...
          +kt*beta*f(i,j+1) ...
          +kb*beta*f(i,j-1) ...
          +Dxs*g(i,j))/(kl+kr+beta*(kb+kt)) ...
          - f(i,j);
%   res = (f(i+1,j)+f(i-1,j)+beta*(f(i,j+1)+f(i,j-1)) ...
%         +Dxs*g(i,j))/beta1 - f(i,j);
   f(i,j) = f(i,j) + relax*res;
   cor = abs(f(i,j)-fsave(i,j));
   if(abs(cor) > cormax)
     cormax = cor;
   end
  end
 end

 cormax

%------------------------------------------------
% wrap: implement the periodic boundary condition
%------------------------------------------------

 for j=1:Ny+1
   f(1,j)    = f(Nx+1,j);
   f(Nx+2,j) = f(2,j);
  end

 if(cormax<tol) 
   Iflag=1;
   break
 end

%---
end % of iterations
%---

%-----
% done
%-----

return
