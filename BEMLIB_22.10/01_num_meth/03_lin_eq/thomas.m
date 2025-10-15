function x = thomas (n,a,b,c,s)

%-----------------------------------------
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%==================================================
% Thomas algorithm for a tridiagonal system
%
% n:     system size
% a,b,c: diagonal, superdiagonal,
%        and subdiagonal elements
% s:     right-hand side
%==================================================

%------------------------------
% reduction to upper bidiagonal
%------------------------------

d(1) = b(1)/a(1);
y(1) = s(1)/a(1);

for i=1:n-1
  i1 = i+1;
   den   = a(i1)-c(i1)*d(i);
   d(i1) = b(i1)/den;
   y(i1) = (s(i1)-c(i1)*y(i))/den;
end

%------------------
% back substitution
%------------------

x(n) = y(n);

for i=n-1:-1:1
  x(i)= y(i)-d(i)*x(i+1);
end

%-----
% done
%-----

return;
