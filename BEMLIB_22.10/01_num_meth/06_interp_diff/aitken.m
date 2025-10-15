function [yint,A] = aitken(N,x,f,xint)

%----
% Aitken interpolation
%----

for i=1:N+1
 u(i) = f(i);
 A(i,1) = u(i);
end

for m=2:N+1
  for i=m:N+1
   u(i) = ((xint-x(i))*u(m-1)-(xint-x(m-1))*u(i))/(x(m-1)-x(i));
   A(i,m) = u(i);
  end
end

yint = u(N+1);

%---
% done
%---

return
