function dydt = blas_fnc(t,y)

%---
% phase-space velocity of the extended Blasius system
%---

 dydt = [y(2);
         y(3);
        -0.5*y(1)*y(3);
         y(5);
         y(6);
        -0.5*(y(4)*y(3)+y(6)*y(1)) ];

%---
% done
%---

 return
