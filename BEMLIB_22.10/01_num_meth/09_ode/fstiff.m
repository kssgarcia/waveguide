function dydt = fstiff(t, y)

%---
% phase-space velocity of a stiff equation
%---

dydt= - 100*(y-sin(t));

return
