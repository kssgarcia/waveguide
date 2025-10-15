function [a,b,c] = splc_clm_true (N,x,f,slope1)

%--------------------------------------
%  Cubic-spline fit of prescribed data
%  with clamped/true end-conditions
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

  at(1) = 2.0D0*(h(1)+h(2))-0.50D0*h(1);
  bt(1) = h(2);

  for i=2:N-1
   at(i) = 2.0D0*(h(i)+h(i+1));
   bt(i) = h(i+1);
   ct(i) = h(i);
  end

  at(N-1)=2.0*(h(N-1)+h(N))+h(N)*(h(N-1)+h(N))/h(N-1);
  ct(N-1)=(h(N-1)^2-h(N)^2)/h(N-1);

  for i=1:N-1
    rhs(i)  = 3.0D0*( (f(i+2)-f(i+1))/h(i+1) ...
                     -(f(i+1)-f(i) )/h(i)  );
  end

  rhs(1) = rhs(1)-1.5*( (f(2)  -f(1))/h(1) - slope1);

%--------------------
% solve N-1 equations
%--------------------

  sln = thomas (N-1,at,bt,ct,rhs);

  for i=1:N-1
   b(i+1) = sln(i);
  end

  b(1) = -0.5D0*b(2) + 1.5D0*((f(2)-f(1))/h(1)-slope1)/h(1);
  b(N+1)=(h(N-1)+h(N))*b(N)/h(N-1)-h(N)*b(N-1)/h(N-1);

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
