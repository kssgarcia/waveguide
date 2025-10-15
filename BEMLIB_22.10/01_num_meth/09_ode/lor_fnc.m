function dydt = lor_fnc(t, y)

%---
% phase-space velocity for the Lorenz system
%---

 k = 10.0;
 b = 8.0/3.0;
 r = 28.0;

 dydt = [-k*(y(1)-y(2));
         -y(1)*y(3)+r*y(1)-y(2);
         y(1)*y(2)-b*y(3)];

%---
% done
%---

return
