 function [x,f,Iflag] = newton2 ...
     ...
    (menu ...
    ,Niter ...
    ,eps ...
    ,x ...
    )

%========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%========================================

%------------------------------------------------
% This program accompanies the book:
%          C. Pozrikidis
% "Numerical Computation in Science and Engineering"
%       Oxford University Press
%------------------------------------------------

%---------------------------------------------
%  Newton's method for two nonlinear equations
%
%  SYMBOLS:
%  --------
%
%  eps:	  small interval for computing the Jacobian
%   	  by numerical differentiation
%  Dx: 	  correction vector
%  tol:	  accuracy
%  Iflag: will set equal to 1 if something is wrong
%--------------------------------------------------

   tol = 0.0000001;
   relax = 1.0;

%-----------
% initialize
%-----------

   Iflag = 1;

%---------------------
% start the iterations
%---------------------

   for Iter=1:Niter

    f = newton2_fun (menu,x);

%---------------------
% compute the Jacobian
% by numerical differentiation
%---------------------

   for j=1:2
     x(j) = x(j)+eps;      % perturb
     f1 = newton2_fun(menu,x);
     x(j) = x(j)-eps;      % reset
       for i=1:2
        Jac(i,j) = (f1(i)-f(i))/eps;
       end
   end

%---
% solve the equation: Jac . Dx = - f
% for the correction vector Dx
% by Cramer's rule
%---

   b1  = -f(1);
   b2  = -f(2);
   det = Jac(1,1)*Jac(2,2)-Jac(1,2)*Jac(2,1);
   dx(1) = (b1*Jac(2,2)-Jac(1,2)*b2)/det;
   dx(2) = (b2*Jac(1,1)-Jac(2,1)*b1)/det;

%--------
% correct
%--------

   x(1) = x(1)+relax*dx(1);
   x(2) = x(2)+relax*dx(2);

%-------
% escape
%-------

  iescape = 1;

  if(abs(dx(1)) > tol) iescape = 0; end
  if(abs(dx(2)) > tol) iescape = 0; end

  if(iescape==1)
    Iflag = 0;
    f = newton2_fun(menu,x);
    return
  end

%----
 end  % of iterations
%----

%-----
% done
%-----

 return
