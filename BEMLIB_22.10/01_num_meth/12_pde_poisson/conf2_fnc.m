function f= conf_map(alpha,k,xi,eta,x,y)

%-------------
% evaluation of the conformal
% mapping functions
%-------------

  cf = 2.0*alpha/k;
  den = cosh(k*y)-cos(k*x);
  f(1) = xi  - x - cf*sin (k*x)/den;
  f(2) = eta - y + cf*sinh(k*y)/den;

%-----
% done
%-----

return
