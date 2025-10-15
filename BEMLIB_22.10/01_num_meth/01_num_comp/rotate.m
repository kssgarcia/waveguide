function rotate

%============================
% Rotation of a point (x) about an axis
% subtended between the origin and the point (X)
%============================

close all
clear all

%---
% define the point X
%---

X = 1.0; Y = 1.0; Z = 1.0;

%---
% direction cosines
%---

norm = sqrt(X^2+Y^2+Z^2);

a = X/norm;
b = Y/norm;
c = Z/norm;

%---
% ancillary
%---

ats = 1.0-a^2; 
bts = 1.0-b^2; 
cts = 1.0-c^2; 

%---
% launch the graphics
%---

figure(1)
plot3([0 X],[0,Y],[0,Z],'r')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
box on
hold on
view(202,14)

%---
% initial point position
%---

x = 0.5;
y = 0.7;
z = 0.5;
xp = [x; y; z];

%---
% number of steps and rotation angle
% at each step
%---

nper=3;
Nsteps = nper*64;
Dphi = nper*2.0*pi/Nsteps;

%---
% differential projection matrix
%---

cs = cos(Dphi);
sn = sin(Dphi);
P  = projection (a,b,c,ats,bts,cts,cs,sn);

%---
% animation
%---

phi = 0.0;

%---------
for rep=1:Nsteps
%---------

figure(1)
hold on

xp = P*xp;
xp = xp-0.002*[1;1;1];

x = xp(1); y = xp(2); z = xp(3);
plot3(x,y,z,'o')

phi = phi+Dphi;
pause(0.1)
hold off

%---
end
%---

return

%==========

function P = projection (a,b,c,ats,bts,cts,cs,sn)

%---
% projection matrix for rotation
%---

csc = 1.0-cs;
P(1,1) = a^2+ats*cs; P(1,2) = a*b*csc;    P(1,3) = a*c*csc;
                     P(2,2) = b^2+bts*cs; P(2,3) = b*c*csc;
                                          P(3,3) = c^2+cts*cs; 
P(2,1) = P(1,2); P(3,1) = P(1,3);  P(3,2) = P(2,3); 

P(1,2) = P(1,2)-c*sn; P(1,3) = P(1,3)+b*sn;
P(2,1) = P(2,1)+c*sn; P(2,3) = P(2,3)-a*sn;
P(3,1) = P(3,1)-b*sn; P(3,2) = P(3,2)+a*sn;

%---
% done
%---

return
