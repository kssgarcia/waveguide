function [a,b,c] = splc_clm (N,x,f,slope1,slope2)

%-----------------------------------------
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%------------------------------------------------
% This program accompanies the book:
%           C. Pozrikidis
% ``Numerical Computation in Science and Engineering''
%      Oxford University Press, 2008
%------------------------------------------------

%--------------------------------------
%  Cubic-spline fit of prescribed data
%  with clamped-end conditions
%
%  SYMBOLS
%  -------
%
%  N  .... number of intervals
%  x .... x coord. of prescribed data
%  f .... y coord. of prescribed data
%
%  a .... polynomial coefficient related to 3rd derivative
%  b .... polynomial coefficient related to 2nd derivative
%  c .... polynomial coefficient related to 1st derivative
%  h .... interval between prescribed data
%
%  a.... diagonal of tridiagonal matrix
%  b.... superdiagonal of tridiagonal matrix
%  c.... subdiagonal of tridiagonal matrix
%
%--------------------------------------

%-----------------------
% compute intervals h(i)
%-----------------------

  for i=1:N
   h(i) = x(i+1)-x(i);
  end

%--------------------------------
% Generate the tridiagonal matrix
% for N-1 equations
%--------------------------------

  a(1) = 2.0D0*(h(1)+h(2))-0.50D0*h(1);
  b(1) = h(2);

  for i=2:N-1
   a(i) = 2.0D0*(h(i)+h(i+1));
   b(i) = h(i+1);
   c(i) = h(i);
  end

  a(N-1) = 2.0D0*(h(N-1)+h(N))-0.50D0*h(N);
  c(N-1) = h(N-1);

  for i=1:N-1
    rhs(i)  = 3.0D0*( (f(i+2)-f(i+1))/h(i+1) ...
                     -(f(i+1)-f(i) )/h(i)  );
  end

 rhs(1)   = rhs(1)   - 1.5*( (f(2)  -f(1))/h(1) - slope1);
 rhs(N-1) = rhs(N-1) + 1.5*( (f(N+1)-f(N))/h(N) - slope2);

%--------------------
% solve N-1 equations
%--------------------

  sln = thomas (N-1,a,b,c,rhs);

  for i=1:N-1
   b(i+1) = sln(i);
  end

  b(1)   = -0.5D0*b(2) + 1.5D0*( (f(2) -f(1))/h(1) - slope1)/h(1);
  b(N+1) = -0.5D0*b(N) - 1.5D0*( (f(N+1)-f(N))/h(N) - slope2)/h(N);

%---------------------------------
% compute the coefficients a and c
%---------------------------------

   for i=1:N
    a(i) = (b(i+1)-b(i))/(3.0D0*h(i));
    c(i) = (f(i+1)-f(i))/h(i) - h(i)*(b(i+1)+2.0D0*b(i))/3.0D0;
   end

%-----
% done
%-----

  return
