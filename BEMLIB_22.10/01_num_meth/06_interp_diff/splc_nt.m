function [a,b,c] = splc_nt (N,x,y)

%-----------------------------------
% natural cubic spline interpolation
%
% N: number of intervals
%-----------------------------------

%----------
% intervals
%----------

for i=1:N
  h(i)=x(i+1)-x(i);
end

%--------------------------------
% generate the tridiagonal matrix
%--------------------------------

at(1) = 2.0*(h(1)+h(2));
bt(1) = h(2);

for i=2:N-2
  at(i) = 2.0*(h(i)+h(i+1));
  bt(i) = h(i+1);
  ct(i) = h(i);
end

at(N-1) = 2.0*(h(N-1)+h(N));
ct(N-1) = h(N-1);

%----------------
% right-hand side
%----------------

for i=1:N-1
   rhs(i) = 3.0*( (y(i+2)-y(i+1))/h(i+1) ...
                 -(y(i+1)-y(i) )/h(i)  );
end

%-----------------
% solve the system
%-----------------

sol = thomas (N-1,at,bt,ct,rhs);

%-----------------
% recover the bees
%-----------------

b(1)=0.0;

for i=1:N-1
  b(i+1) = sol(i);
end

b(N+1)=0.0;

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
