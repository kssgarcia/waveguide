function [f,iter,Iflag] = pois_gs_nnnn ...
...
   (Nx,Ny,Dx,Dy,g,itermax,tol,relax ...
...
   ,qleft,qright,qbot,qtop,f,Ishift)

%=============================================
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
% The solution is found by point 
% Gauss-Siedel iterations
%
% The source at the interior and boundary nodes
% is used
%=============================================

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

method = 1;

%------------------------
% Gauss-Siedel iterations
%------------------------

for iter=1:itermax

%------------------------
% update nodes row-by-row
%------------------------

 fsv = f;

%---
% interior nodes
%---

 for j=2:Ny
  for i=2:Nx

   if(method==0)
    res = (f(i+1,j)+f(i-1,j)+beta*(f(i,j+1)+f(i,j-1)) ...
          +Dxs*g(i,j))/beta1 - f(i,j);
    f(i,j) = f(i,j) + relax*res;
   elseif(method==1)
    res = (fsv(i+1,j)+fsv(i-1,j)+beta*(fsv(i,j+1)+fsv(i,j-1)) ...
          +Dxs*g(i,j))/beta1 - fsv(i,j);
    f(i,j) = fsv(i,j) + relax*res;
   end

  end
 end

%-----
% left
%-----

 i=1;
 for j=2:Ny

   if(method==0)
   res = (2*fsv(i+1,j)-Dx2*qleft+beta*(fsv(i,j+1)+fsv(i,j-1)) ...
                +Dxs*g(i,j))/beta1 - fsv(i,j);
    f(i,j) = fsv(i,j) + relax*res;
   elseif(method==1)
    res = (2*f(i+1,j)-Dx2*qleft+beta*(f(i,j+1)+f(i,j-1)) ...
                 +Dxs*g(i,j))/beta1 - f(i,j);
     f(i,j) = f(i,j) + relax*res;
   end

 end

% corner points:

 j=1;
 if(method==0)
  res = (2*fsv(i+1,j)-Dx2*qleft ...
        +beta*(fsv(i,j+1)+fsv(i,j+1)-Dy2*qbot) ...
               +Dxs*g(i,j))/beta1 - fsv(i,j);
  f(i,j) = fsv(i,j) + relax*res;
 elseif(method==1)
  res = (2*f(i+1,j)-Dx2*qleft+beta*(f(i,j+1)+f(i,j+1)-Dy2*qbot) ...
               +Dxs*g(i,j))/beta1 - f(i,j);
  f(i,j) = f(i,j) + relax*res;
 end

 j=Ny+1;
 if(method==0)
  res = (2*fsv(i+1,j)-Dx2*qleft ...
       +beta*(fsv(i,j-1)+fsv(i,j-1)-Dy2*qtop) ...
               +Dxs*g(i,j))/beta1 - fsv(i,j);
  f(i,j) = fsv(i,j) + relax*res;
 elseif(method==1)
  res = (2*f(i+1,j)-Dx2*qleft+beta*(f(i,j-1)+f(i,j-1)-Dy2*qtop) ...
               +Dxs*g(i,j))/beta1 - f(i,j);
  f(i,j) = f(i,j) + relax*res;
 end
%---

%-----
% right
%-----

 i=Nx+1;
 for j=2:Ny
  if(method==0)
   res = (2*fsv(i-1,j)-Dx2*qright ...
      +beta*(fsv(i,j+1)+fsv(i,j-1)) ...
                +Dxs*g(i,j))/beta1 - fsv(i,j);
   f(i,j) = fsv(i,j) + relax*res;
  elseif(method==1)
   f(i,j) = fsv(i,j) + relax*res;
   res = (2*f(i-1,j)-Dx2*qright+beta*(f(i,j+1)+f(i,j-1)) ...
                +Dxs*g(i,j))/beta1 - f(i,j);
   f(i,j) = f(i,j) + relax*res;
  end
 end

% corner points:

 j=1;
 if(method==0)
  res = (2*fsv(i-1,j)-Dx2*qright ...
      +beta*(fsv(i,j+1)+fsv(i,j+1)-Dy2*qbot) ...
               +Dxs*g(i,j))/beta1 - fsv(i,j);
  f(i,j) = fsv(i,j) + relax*res;
 elseif(method==1)
  res = (2*f(i-1,j)-Dx2*qright+beta*(f(i,j+1)+f(i,j+1)-Dy2*qbot) ...
               +Dxs*g(i,j))/beta1 - f(i,j);
  f(i,j) = f(i,j) + relax*res;
 end

 j=Ny+1;
 if(method==0)
  res = (2*fsv(i-1,j)-Dx2*qright ...
     +beta*(fsv(i,j-1)+fsv(i,j-1)-Dy2*qtop) ...
               +Dxs*g(i,j))/beta1 - fsv(i,j);
  f(i,j) = fsv(i,j) + relax*res;
 elseif(method==1)
  res = (2*f(i-1,j)-Dx2*qright+beta*(f(i,j-1)+f(i,j-1)-Dy2*qtop) ...
               +Dxs*g(i,j))/beta1 - f(i,j);
  f(i,j) = f(i,j) + relax*res;
 end

%-----
% bottom
%-----

 j=1;
 for i=2:Nx
   if(method==0)
   res = (fsv(i+1,j)+fsv(i-1,j) ...
    +beta*(2*fsv(i,j+1)-Dy2*qbot) ...
                +Dxs*g(i,j))/beta1 - fsv(i,j);
   f(i,j) = fsv(i,j) + relax*res;
   elseif(method==1)
   res = (f(i+1,j)+f(i-1,j)+beta*(2*f(i,j+1)-Dy2*qbot) ...
                +Dxs*g(i,j))/beta1 - f(i,j);
   f(i,j) = f(i,j) + relax*res;
   end
 end

%-----
% top
%-----

 j=Ny+1;
 for i=2:Nx
   if(method==0)
   res = (fsv(i+1,j)+fsv(i-1,j) ...
     +beta*(2*fsv(i,j-1)-Dy2*qtop) ...
                +Dxs*g(i,j))/beta1 - fsv(i,j);
   f(i,j) = fsv(i,j) + relax*res;
   elseif(method==1)
    res = (f(i+1,j)+f(i-1,j)+beta*(2*f(i,j-1)-Dy2*qtop) ...
                 +Dxs*g(i,j))/beta1 - f(i,j);
    f(i,j) = f(i,j) + relax*res;
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

%-----
% done
%-----

return
