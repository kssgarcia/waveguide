function [a,b,c] = splc_pr (N,x,f)

%======================================================
% This program accompanies the book:
%           C. Pozrikidis
% ``Numerical Computation in Science and Engineering''
%      Oxford University Press
%======================================================

%--------------------------------------
%  Cubic spline interpolation of prescribed data
%  with periodic boundary conditions
%
%  ith cubic: f = a x^3 + b x^2 + c x + f_i
%
%  SYMBOLS
%  -------
%
%  N ... number of intervals
%  x ... x coordinates of prescribed data
%  f ... function values of prescribed data
%
%  c ... polynomial coefficient related to 1st derivative
%  b ... polynomial coefficient related to 2nd derivative
%  a ... polynomial coefficient related to 3rd derivative
%  h ... intervals between prescribed data
%--------------------------------------

%----------------------
% compute the intervals
%----------------------

 for i=1:N
   h(i) = x(i+1)-x(i);
 end

%--------------------------------
% Generate the bordered tridiagonal matrix:
%
% | M11 M12  0    0  . . .  0            M1N   |
% | M21 M22  M23  0  . . .  0            0     |
% | 0   M32  M33 M34 . . .  0            0     |
% | 0   0 ...               0            0     |
% | 0   0 ... ...         M(N-1)(N-1)  M(N-1)N |
% | MN1 0 ...  0           MN(N-1)      MNN    |
%
%--------------------------------

   M = zeros(N,N);

   M(1,1) = 2.0D0*(h(1)+h(2));
   M(1,2) = h(2);
   M(1,N) = h(1);

   for i=2:N-1
     M(i,i)   = 2.0D0*(h(i)+h(i+1));
     M(i,i+1) = h(i+1);
     M(i,i-1) = h(i);
    end

   M(N,N)   = 2.0D0*(h(N)+h(1));
   M(N,1)   = h(1);
   M(N,N-1) = h(N);

%-----
% generate the right-hand side
%-----

  for i=1:N-1
     rhs(i) = 3.0D0*( (f(i+2)-f(i+1))/h(i+1) - (f(i+1)-f(i) )/h(i) );
  end
  rhs(N) = 3.0D0*( (f(2)-f(1))/h(1) - (f(N+1)-f(N))/h(N) );

%-----
% solve the system
%-----

 sln = rhs/M';

%-------------------------
% extract the coefficients
%-------------------------

 for i=1:N
   b(i+1) = sln(i);
 end

 b(1) = b(N+1);

%---------------------------------
% Compute the coefficients a and c
%---------------------------------

for i=1:N
  a(i) = (b(i+1)-b(i))/(3.0D0*h(i));
  c(i) = (f(i+1)-f(i))/h(i) - h(i)*(b(i+1)+2.0D0*b(i))/3.0D0;
end

%-----
% done
%-----

return
