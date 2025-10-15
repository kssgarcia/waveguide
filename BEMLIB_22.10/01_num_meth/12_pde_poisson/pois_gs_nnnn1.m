function [f,iter,Iflag] = fpois_gs_nnnn1 ...
...
   (Nx,Ny,Dx,Dy,g,itermax,tol,relax ...
   ,qleft,qright,qbot,qtop,f,Ishift)

%===============================================
% Solution of Poisson's equation
% in a rectangular domain
% with the uniform Neumann boundary condition
% along the four sides:
%
% bottom:  df/dy =  qbot
% top:     df/dy = -qtop
% left:    df/dx =  qleft
% right:   df/dx = -qright
%
% The solution is found by
% point Gauss--Siedel iterations
%
% The source at the interior nodes is only used
%==============================================

%--------
% prepare
%--------

Dx2 = 2.0*Dx;
Dy2 = 2.0*Dy;

Dxs = Dx*Dx;
Dys = Dy*Dy;

beta = Dxs/Dys;
beta1 = 2.0*(1.0+beta);

Iflag = 0;  % convergence flag, 1 indicates convergence

%------------------------
% Gauss-Siedel iterations
%------------------------

for iter=1:itermax

%------------------------
% update nodes row-by-row
%------------------------

 fsv = f;

%---
% left and right boundaries
%---

 for j=2:Ny
  f(1,j) = (4.0*f(2,j)-f(3,j)-Dx2*qleft)/3.0;
  f(Nx+1,j) = (4.0*f(Nx,j)-f(Nx-1,j)-Dx2*qright)/3.0;
 end

%---
% bottom and top boundaries
%---

 for i=2:Nx
  f(i,1) = (4.0*f(i,2)-f(i,3)-Dy2*qbot)/3.0;
  f(i,Ny+1) = (4.0*f(i,Ny)-f(i,Ny-1)-Dy2*qtop)/3.0;
 end

%---
% interior nodes
%---

method = 1;

  for j=2:Ny
   for i=2:Nx
    if(method==0)
     res = (fsv(i+1,j)+fsv(i-1,j)+beta*(fsv(i,j+1)+fsv(i,j-1)) ...
           +Dxs*g(i,j))/beta1 - fsv(i,j);
     f(i,j) = fsv(i,j) + relax*res;
    elseif(method==1)
     res = (f(i+1,j)+f(i-1,j)+beta*(f(i,j+1)+f(i,j-1)) ...
           +Dxs*g(i,j))/beta1 - f(i,j);
     f(i,j) = f(i,j) + relax*res;
    end
   end
  end

%------
% shift
%------

  if(Ishift==1)

  shift = f(Nx/2,Ny/2);

  for i=1:Nx+1
   for j=1:Ny+1
    f(i,j)=f(i,j)-shift;
   end
  end

  end

%------
% correction
%------

  cormax=0;

  for i=1:Nx+1
   for j=1:Ny+1
   cor=abs(f(i,j)-fsv(i,j));
   if(abs(cor) > cormax)
     cormax = cor;
   end
   end
  end

%-----
% stopping check
%-----

  cormax

  if(cormax<tol) 
   Iflag=1;
   break
  end

%---
end % of iterations
%---

%---
% corner nodes
%---

 f(1,1) = 0.5*(f(1,2)+f(2,1));
 f(Nx+1,1) = 0.5*(f(Nx,1)+f(Nx+1,2));
 f(1,Ny+1) = 0.5*(f(1,Ny)+f(2,Ny+1));
 f(Nx+1,Ny+1) = 0.5*(f(Nx+1,Ny)+f(Nx,Ny+1));

%-----
% done
%-----

return
