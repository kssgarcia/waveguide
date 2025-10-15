function [f, iter, Iflag] = fpois_gs_dpr ...
...
   (Nx,Ny,Dx,Dy,g,itermax,tol,relax,fbot,ftop)

%------------------------------------------
% Solution of Poisson's equation
% in a rectangular domain
% with a uniform Dirichlet boundary condition
% at the top and bottom boundaries
% and periodic boundary condition
% at the left and right boundaries
%
% The solution is found by point Gauss--Siedel
% iterations
%------------------------------------------

%--------
% prepare
%--------

Dxs = Dx*Dx;
Dys = Dy*Dy;
beta = Dxs/Dys;
beta1 = 2.0*(beta+1.0);

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
   res = (f(i+1,j)+f(i-1,j)+beta*(f(i,j+1)+f(i,j-1)) ...
         +Dxs*g(i,j))/beta1 - f(i,j);
   f(i,j) = f(i,j) + relax*res;
   cor=abs(f(i,j)-fsave(i,j));
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
