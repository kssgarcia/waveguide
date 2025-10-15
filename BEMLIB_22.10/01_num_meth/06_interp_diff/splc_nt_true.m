function [a,b,c] = splc_nt_true (N,x,y)

%===============================================
% FDLIB BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%===============================================

%--------------------------------
% cubic spline interpolation with a
% natural-end condition and a true-end condition
%
% N+1 points are provided
%--------------------------------

%----------
% intervals
%----------

for i=1:N
  h(i)=x(i+1)-x(i);
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

at(N-1) = 2.0*(h(N-1)+h(N))+h(N)*(h(N-1)+h(N))/h(N-1);
ct(N-1)=(h(N-1)^2-h(N)^2)/h(N-1);

%----------------
% right-hand side
%----------------

for i=1:N-1
 rhs(i) = 3.0*( (y(i+2)-y(i+1))/h(i+1) ...
               -(y(i+1)-y(i) )/h(i) );
end

%------------------------------
% solve the triadiagonal system
%------------------------------

sol = thomas (N-1,at,bt,ct,rhs);

%-----------------
% recover the bees
%-----------------

for i=1:N-1
  b(i+1) = sol(i);
end

b(1) = 0.0;
b(N+1) = (h(N-1)+h(N))*b(N)/h(N-1)-h(N)*b(N-1)/h(N-1);

%---------------------------------
% recover the coefficients a and c
%---------------------------------

for i=1:N
  a(i) = (b(i+1)-b(i))/(3.0D0*h(i));
  c(i) = (y(i+1)-y(i))/h(i)-h(i)*(b(i+1)+2.0*b(i))/3.0;
end

%---
% done
%---

return
