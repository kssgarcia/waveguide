function f = newton1_2_fun(menu,x)

global rlp
global mexp alpha beta

if(menu==110)
  rj0 = besselj(0,x);
  rj1 = besselj(1,x);
  f = 1.0 + (rj1/rj0-2.0/x)*rj1/rj0;
  f = f-2.0*rlp;
elseif(menu==111)
  f = abs(x+alpha)^mexp - abs(x-alpha)^mexp-beta;
elseif(menu==200)
  f = 2.0*sin(pi*x/2)^2 - 0.5*cosh(0.5*x)^3;
end

%---
% done
%---

return
