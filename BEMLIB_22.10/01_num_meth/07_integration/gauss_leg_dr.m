function integral = gauss_leg_dr(a,b,NQ)

%====================================
% Driver for Gauss-Legendre quadrature
% will integrate from x=a to b
%
% NQ: number of base points
%====================================

%---
% prepare
%---

[Z,W] = gauss_leg(NQ);

xm = 0.5*(b+a);
xd = 0.5*(b-a);

%---
% launch
%---

integral = 0.0;

for q=1:NQ
 x = xm+xd*Z(q);
 fnc = x;
 integral = integral + fnc*W(q);
end

integral = integral*xd;

%---
% done
%---

return
