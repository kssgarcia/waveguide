function [a,b,c] = splc_true(N,x,y)

%================================
% FDLIB
%
% true cubic spline interpolation
% if a function y=f(x)
%
% if the data fall on a cubic,
% the interpolation is exact
%
% N: number of intervals
%================================

%------------------------------------------------
% This program accompanies the book:
%           C. Pozrikidis
% ``Numerical Computation in Science and Engineering''
%      Oxford University Press
%------------------------------------------------

%----------
% intervals
%----------

for i=1:N
  h(i) = x(i+1)-x(i);
end

%------------------------------
% generate a tridiagonal matrix
%------------------------------

at(1) = 2.0*(h(1)+h(2));
bt(1) = h(2);

for i=2:N-2
  at(i) = 2.0*(h(i)+h(i+1));
  bt(i) = h(i+1);
  ct(i) = h(i);
end

at(1) = at(1)+h(1)*(h(1)+h(2))/h(2);
bt(1) = (h(2)^2-h(1)^2)/h(2);

at(N-1) = 2.0*(h(N-1)+h(N))+h(N)*(h(N-1)+h(N))/h(N-1);
ct(N-1) = (h(N-1)*h(N-1)-h(N)*h(N))/h(N-1);

%----------------
% right-hand side
%----------------

for i=1:N-1
 rhs(i) = 3.0*( (y(i+2)-y(i+1))/h(i+1) ...
               -(y(i+1)-y(i) )/h(i) );
end

%-----------------
% solve the system
%-----------------

sol = thomas(N-1,at,bt,ct,rhs);

%-----------------
% recover the bees
%-----------------

for i=1:N-1
  b(i+1) = sol(i);
end

b(1) = (h(1)+h(2))*b(2)/h(2) -h(1)*b(3)/h(2);
b(N+1) = (h(N-1)+h(N))*b(N)/h(N-1)-h(N)*b(N-1)/h(N-1);

%---------------------
% coefficients a and c
%---------------------

for i=1:N
  a(i) = (b(i+1)-b(i))/(3.0D0*h(i));
  c(i) = (y(i+1)-y(i))/h(i)-h(i)*(b(i+1)+2.0*b(i))/3.0;
end

%---
% done
%---

return
