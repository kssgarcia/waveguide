clear all
close all

Ndiv=16;

%---
% animation of the rotation of a disk
%---

box
view(-27, 46)

%---
% disk in the yz plane
%---

%[YREF,ZREF,XREF]=sphere(Ndiv);
%XREF=0.1*XREF;

for i=1:Ndiv+1
 theta = 2*pi*(i-1.0)/Ndiv;
 XREF(i)=0.0;
 YREF(i)= cos(theta);
 ZREF(i)= sin(theta);
end

%---
% first point of the director is fixed
%---

xplt(1)=0.0;
yplt(1)=0.0;
zplt(1)=0.0;

Dt = 0.02;

%---
for iloop=1:20000
%---

time=(iloop-1.0)*Dt;

dx= 1.0+2.0*cos(time);
dy= 1.0+sin(time);
dz= 1.0+5*sin(time);

%---
% unit director
%---

dm=sqrt(dx^2+dy^2+dz^2);
dx=dx/dm;
dy=dy/dm;
dz=dz/dm;

xplt(2)=dx;
yplt(2)=dy;
zplt(2)=dz;

%---
% deduce the meridional and azimuthal angles
%---

ds=sqrt(dy^2+dz^2);
dr=sqrt(dx^2+dy^2+dz^2);

csph = dy/ds;
snph = dz/ds;
phi  = acos(dy/dr);
if(dz<0)
  phi=2*pi-phi;
end

csth =dx/dr;
snth = sqrt(1.0-csth^2);
ROTZ=eye(3);
ROTZ(1,1)=csth;
ROTZ(2,2)=csth;
ROTZ(1,2)=-snth;
ROTZ(2,1)= snth;

ROTX=eye(3);
ROTX(2,2)= csph;
ROTX(3,3)= csph;
ROTX(2,3)=-snph;
ROTX(3,2)= snph;

%for i=1:Ndiv+1
% for j=1:Ndiv+1
% v(1)=XREF(i,j);
% v(2)=YREF(i,j);
% v(3)=ZREF(i,j);
% u=ROTX*ROTZ*v';
% X(i,j)=u(1);
% Y(i,j)=u(2);
% Z(i,j)=u(3);
% end
%end

for i=1:Ndiv+1
  v(1)=XREF(i);
  v(2)=YREF(i);
  v(3)=ZREF(i);
  u=ROTX*ROTZ*v';
  X(i)=u(1);
  Y(i)=u(2);
  Z(i)=u(3);
end

if(iloop==1)
% Handle1 = surf(X,Y,Z);
  Handle1 = patch(X,Y,Z);
 set(Handle1,'EraseMode','xor')
 hold on
 Handle2 =plot3(xplt,yplt,zplt);
% set(Handle2,'EraseMode','xor')
 axis([-2 2 -2 2 -2 2])
 %axis square;
 xlabel('x')
 ylabel('y')
 zlabel('z')
end

 set(Handle1,'XData',X,'YData',Y,'ZData',Z);
 set(Handle2,'XData',xplt,'YData',yplt,'ZData',zplt);
 drawnow

%---
end
%---
