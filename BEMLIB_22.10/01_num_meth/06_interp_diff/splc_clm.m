function [a,b,c] = splc_clm (N,x,f,slopeL,slopeR)

%-----------------------------------------
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%----------------------------------------

%------------------------------------------------
% This program accompanies the book:
%           C. Pozrikidis
% ``Numerical Computation in Science and Engineering''
%      Oxford University Press
%------------------------------------------------

%--------------------------------------
%  Cubic-spline fit of prescribed data
%  f(x) with clamped-end conditions
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
%--------------------------------------

%-----------------------
% compute intervals h(i)
%-----------------------

  for i=1:N
   h(i) = x(i+1)-x(i);
  end

%--------------------------------
% generate a tridiagonal matrix
% for N-1 equations
%--------------------------------

  at(1) = 2.0D0*(h(1)+h(2))-0.50D0*h(1);
  bt(1) = h(2);

  for i=2:N-1
   at(i) = 2.0D0*(h(i)+h(i+1));
   bt(i) = h(i+1);
   ct(i) = h(i);
  end

  at(N-1) = 2.0D0*(h(N-1)+h(N))-0.50D0*h(N);
  ct(N-1) = h(N-1);

  for i=1:N-1
    rhs(i)  = 3.0D0*( (f(i+2)-f(i+1))/h(i+1) ...
                     -(f(i+1)-f(i) )/h(i)  );
  end

 rhs(1)   = rhs(1)   - 1.5*( (f(2)  -f(1))/h(1)-slopeL);
 rhs(N-1) = rhs(N-1) + 1.5*( (f(N+1)-f(N))/h(N)-slopeR);

%--------------------
% solve N-1 equations
%--------------------

  sln = thomas (N-1,at,bt,ct,rhs);

%-----------------
% recover the bees
%-----------------

  for i=1:N-1
   b(i+1) = sln(i);
  end

  b(1)  =-0.5D0*b(2)+1.5*((f(2)  -f(1))/h(1)-slopeL)/h(1);
  b(N+1)=-0.5D0*b(N)-1.5*((f(N+1)-f(N))/h(N)-slopeR)/h(N);

%---------------------------------
% compute the coefficients a and c
%---------------------------------

 for i=1:N
   a(i) = (b(i+1)-b(i))/(3.0D0*h(i));
   c(i) = (f(i+1)-f(i))/h(i)-h(i)*(b(i+1)+2.0D0*b(i))/3.0D0;
 end

%-----
% done
%-----

return
