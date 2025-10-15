close all
clear all

%------------------------------------------------
% This program accompanies the book:
%
%          C. Pozrikidis
% ``Numerical Computation in Science and Engineering''
%         Oxford University Press, 1998
%------------------------------------------------
 
%--------------------------------------
% anination of a 3D line interpolated by cubic splines
%
% LEGEND:
% ------
%
% N:        Number of intervals
% xp(i):    x value of ith point
% fp(i):    f value of ith point
% a, b, c:  cubic spline coefficients
%
% kpl:  number of intervals for plotting the spline
%
%===============================================

%---
% node coordinates
%---

N = 10;

fpx = [0.0  0.1  0.2  0.3  0.4  0.5  0.4  0.3  0.2  0.1  ];
fpy = [0.01 0.12 0.43 0.34 0.45 0.54 0.13 0.32 0.21 0.11 ];
fpz = [0.02 0.14 0.26 0.38 0.46 0.54 0.42 0.30 0.22 0.14 ];

%---
% initialize the node velocity
%---

for i=1:N
 velx(i) = 1.0;
 vely(i) = 0.9;
 velz(i) = 0.8;
end

%---
% animation loop
%---

for repeat=1:9000

%---
% points move randomly
%---

Ido = 1;
Ido = 0;

 if(Ido==1)
  Dt = 0.01;
  for i=1:N
   fpx(i) = fpx(i)+Dt*(rand(1)-0.5);
   fpy(i) = fpy(i)+Dt*(rand(1)-0.5);
   fpz(i) = fpz(i)+Dt*(rand(1)-0.5);
  end
 end

%---
% points bounce in a box
%---

Ido = 0;
Ido = 1;

 if(Ido==1)
  Dt = 0.01;
  win =0.5;
  for i=1:N
   fpx(i) = fpx(i)+velx(i)*Dt;
   fpy(i) = fpy(i)+vely(i)*Dt;
   fpz(i) = fpz(i)+velz(i)*Dt;
   if(fpx(i)> win) velx(i)=-abs(velx(i)); end
   if(fpx(i)<-win) velx(i)= abs(velx(i)); end
   if(fpy(i)> win) vely(i)=-abs(vely(i)); end
   if(fpy(i)<-win) vely(i)= abs(vely(i)); end
   if(fpz(i)> win) velz(i)=-abs(velz(i)); end
   if(fpz(i)<-win) velz(i)= abs(velz(i)); end
  end
  axis([-win win -win win -win win])
end

%-----
% wrap
%-----

 fpx(N+1) = fpx(1); 
 fpy(N+1) = fpy(1); 
 fpz(N+1) = fpz(1); 

%====
% compute the polygonal arc length
% to be used as the interpolation variable
%====
 
xp(1) = 0.0;
 
for i=2:N+1
 xp(i) = xp(i-1)+ sqrt((fpx(i)-fpx(i-1))^2+(fpy(i)-fpy(i-1))^2 ...
                      +(fpz(i)-fpz(i-1))^2);
end

%====
% spline coefficients
%====

[ax,bx,cx] = splc_pr(N,xp,fpx);
[ay,by,cy] = splc_pr(N,xp,fpy);
[az,bz,cz] = splc_pr(N,xp,fpz);

%====
% plotting
%====

kpl = 128; % number of plotting intervals

dxint = (xp(N+1)-xp(1))/kpl;
 
%===
% loop over plotting points
%===

xint = xp(1)+0.00001; % initial point
 
for j=1:kpl+1
 
%---
% find the host interval
%---

 Iskip = 0;

 for i=1:N
  prod = (xint-xp(i))*(xint-xp(i+1));
    if(prod < 0.00)
       l = i;
       Iskip = 1;
     end
  if(Iskip==1) break; end
 end
 
%---
% compute the spline
%---

 xd = xint-xp(l);
 
 xx(j) = ( (ax(l)*xd+ bx(l) )*xd+ cx(l) )*xd + fpx(l);
 yy(j) = ( (ay(l)*xd+ by(l) )*xd+ cy(l) )*xd + fpy(l);
 zz(j) = ( (az(l)*xd+ bz(l) )*xd+ cz(l) )*xd + fpz(l);

 xint = xint+dxint;

end

%---
% animation
%---

if(repeat==1) 
 Handle1 = plot3(xx,yy,zz,'-','linewidth',3);
 box
 hold on
 Handle2 = plot3(fpx,fpy,fpz,'ro');
 set(Handle1,'EraseMode','xor')
 set(Handle2,'EraseMode','xor')
 xlabel('x')
 ylabel('y')
 zlabel('z')
end

set(Handle1,'XData',xx,'YData',yy,'ZData',zz);
set(Handle2,'XData',fpx,'YData',fpy,'ZData',fpz)
drawnow 

%--
end % of animation loop
%--
