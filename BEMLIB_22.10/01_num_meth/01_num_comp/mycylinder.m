function mycylinder (x1,y1,z1,x2,y2,z2,rad,Ndiv)

%=====================================
% draw a circular cylinder of radius rad
% whose centerline begins at x1, y1, z1
% and ends at x2, y2, z2
% periphery is defined by Ndiv divisions
%=====================================

%---
% unit vector
%---

vx=x2-x1; vy=y2-y1; vz=z2-z1;
vm=sqrt(vx*vx+vy*vy+vz*vz);
vx=vx/vm; vy=vy/vm; vz=vz/vm;

%---
% first rotate by theta
%---

theta=acos(vx);

%---
% second rotate by phi
%---

vn=sqrt(vy*vy+vz*vz);
phi = acos(vy/vn);
if(vz<0)
 phi = 2*pi-phi;
end

%---
% rotation matrices
%---

cs=cos(theta);
sn=sin(theta);
R1=[cs -sn 0;
    sn  cs 0;
    0 0 1];

cs=cos(phi);
sn=sin(phi);
R2=[1 0 0;
    0 cs -sn;
    0 sn cs];
R=R2*R1;

%---
% strips around the cylinder
%---

dtheta=2*pi/Ndiv;

for i=1:Ndiv+1
 x=0.0;
 theta=(i-1)*dtheta;
 y = rad*cos(theta);
 z = rad*sin(theta);
 xx1=R*[x;y;z];
 xcyl(i)=xx1(1);
 ycyl(i)=xx1(2);
 zcyl(i)=xx1(3);
end

xcyl1=xcyl+x1;
ycyl1=ycyl+y1;
zcyl1=zcyl+z1;
xcyl2=xcyl+x2;
ycyl2=ycyl+y2;
zcyl2=zcyl+z2;

%---
% plot patches
%---

hold on
plot3([x1 x2],[y1 y2],[z1 z2],'r')

for i=1:Ndiv
patch([xcyl1(i), xcyl2(i), xcyl2(i+1), xcyl1(i+1), xcyl1(i) ]...
     ,[ycyl1(i), ycyl2(i), ycyl2(i+1), ycyl1(i+1), ycyl1(i) ]...
     ,[zcyl1(i), zcyl2(i), zcyl2(i+1), zcyl1(i+1), zcyl1(i) ]...
     ,'c');
end

patch(xcyl1,ycyl1,zcyl1,'y')
patch(xcyl2,ycyl2,zcyl2,'y')

xlabel('x')
ylabel('y')
zlabel('z')
box
axis equal

%---
% done
%---

return
