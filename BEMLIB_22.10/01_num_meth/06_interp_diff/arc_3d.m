function [xc,yc,zc,a,chi1,chi3] = arc_3d ...
...
      (x1,x2,x3,y1,y2,y3,z1,z2,z3)

%===========================================
% Compute a circular arc passing through
% three points in space
%
% (xc, yc, zc): arc center
%
% chi1: angle of the first point in the plane of the arc
%       measured with respect to the second point
%
% chi3: angle of the third point in the plane of the arc
%       measured with respect to the second point
%===========================================

M(1,1) = 2.0*(x1-x2);
M(1,2) = 2.0*(y1-y2);
M(1,3) = 2.0*(z1-z2);
M(2,1) = 2.0*(x3-x2);
M(2,2) = 2.0*(y3-y2);
M(2,3) = 2.0*(z3-z2);

tmp = x2^2+y2^2+z2^2;

b(1) = x1^2 + y1^2+z1^2 - tmp;
b(2) = x3^2 + y3^2+z3^2 - tmp;

crx = (y3-y1)*(z2-z1)-(z3-z1)*(y2-y1);
cry = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1);
crz = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);
M(3,1) = crx;
M(3,2) = cry;
M(3,3) = crz;
b(3) = x1*crx + y1*cry + z1*crz

sol = b/M';
xc = sol(1); yc=sol(2); zc=sol(3);

as = (x1-xc)^2+(y1-yc)^2+(z1-zc)^2;
a  = sqrt(as);
prj1 = (x1-xc)*(x2-xc)+(y1-yc)*(y2-yc)+(z1-zc)*(z2-zc);
prj3 = (x3-xc)*(x2-xc)+(y3-yc)*(y2-yc)+(z3-zc)*(z2-zc);

chi1 = - acos(prj1/as);
chi3 =   acos(prj3/as);

%---
% done
%---

return
