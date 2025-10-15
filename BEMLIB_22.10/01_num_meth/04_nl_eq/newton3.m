  function [x,f,Iflag] = newton3 ...
     ...
    (menu ...
    ,Niter ...
    ,eps ...
    ,x ...
    )

%====================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%====================================

%------------------------------------------------
% This program accompanies the book:
%
%          C. Pozrikidis
% "Numerical Computation in Science and Engineering"
%       Oxford University Press
%------------------------------------------------

%---------------------------------------------
%  Newton's method for three nonlinear equations
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

    f = newton3_fun(menu,x);

%---------------------
% compute the Jacobian
% by numerical differentiation
%---------------------

      for j=1:3
        x(j) = x(j)+eps;      % perturb
        f1 = newton3_fun(menu,x);
        x(j) = x(j)-eps;      % reset
        for i=1:3
          Jac(i,j) = (f1(i)-f(i))/eps;
        end
      end

%---
% solve the equation: Jac . Dx = - f
% for the correction vector Dx
% by Cramer's rule
%---

      A11 = Jac(1,1); A12 = Jac(1,2); A13 = Jac(1,3);
      A21 = Jac(2,1); A22 = Jac(2,2); A23 = Jac(2,3);
      A31 = Jac(3,1); A32 = Jac(3,2); A33 = Jac(3,3);

      B1  = -f(1);
      B2  = -f(2);
      B3  = -f(3);

      Det =  A11*( A22*A33-A23*A32 ) ...
           - A12*( A21*A33-A23*A31 ) ...
           + A13*( A21*A32-A22*A31 );

      Det1 =  B1*( A22*A33-A23*A32 ) ...
           - A12*(  B2*A33-A23*B3  ) ...
           + A13*(  B2*A32-A22*B3  );

      Det2 = A11*( B2 *A33-A23*B3  ) ...
           -  B1*( A21*A33-A23*A31 ) ...
           + A13*( A21* B3-B2 *A31 );

      Det3 = A11*( A22* B3-A32* B2 ) ...
           - A12*( A21* B3-A31* B2 ) ...
           +  B1*( A21*A32-A22*A31 );

      dx(1) = Det1/Det;
      dx(2) = Det2/Det;
      dx(3) = Det3/Det;

%--------
% correct
%--------

     x(1) = x(1)+relax*dx(1);
     x(2) = x(2)+relax*dx(2);
     x(3) = x(3)+relax*dx(3);

%-------
% escape
%-------

      Iescape = 1;
      if(abs(dx(1)) > tol) Iescape = 0; end
      if(abs(dx(2)) > tol) Iescape = 0; end
      if(abs(dx(3)) > tol) Iescape = 0; end

      if(Iescape==1)
       Iflag = 0;
       f = newton3_fun(menu,x);
       return
      end

%----
     end  % of iterations
%----

%-----
% done
%-----

 return
