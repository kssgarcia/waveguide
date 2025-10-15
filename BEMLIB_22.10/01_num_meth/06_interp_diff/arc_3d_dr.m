%---
clear all
close all
%---

%===========================================
% Compute a circular arc passing through
% three points in space
%===========================================

x1=0.2; y1=0.2; z1=0.0;
x2=1.0; y2=0.3; z2=0.1;
x3=1.3; y3=0.5; z3=0.5;

[xc,yc,zc,a,chi1,chi3] = arc_3d ...
...
      (x1,x2,x3,y1,y2,y3,z1,z2,z3)

npl=32;
Dchi  = (chi3-chi1)/npl;
as = a*a;

M(1,1) = x1-xc; M(1,2) = y1-yc; M(1,3) = z1-zc;
M(2,1) = x2-xc; M(2,2) = y2-yc; M(2,3) = z2-zc;
M(3,1) = (y3-y1)*(z2-z1)-(z3-z1)*(y2-y1);
M(3,2) = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1);
M(3,3) = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);

for i=1:npl+1

 chi = chi1+(i-1.0)*Dchi;
 b(1) = as*cos(chi-chi1);
 b(2) = as*cos(chi);
 b(3) = 0.0;
 sol = b/M';
 x(i) = sol(1) +xc;
 y(i) = sol(2) +yc;
 z(i) = sol(3) +zc;

end

%---
figure(1)
hold on
plot3(x1,y1,z1,'o')
plot3(x2,y2,z2,'d')
plot3(x3,y3,z3,'s')
plot3(x,y,z)
axis square
box
xlabel('x','fontsize',15)
ylabel('y','fontsize',15)
zlabel('z','fontsize',15)
set(gca,'fontsize',15)
%---

