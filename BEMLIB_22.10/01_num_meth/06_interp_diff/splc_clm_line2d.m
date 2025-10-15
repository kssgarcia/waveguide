%------------------------------------------------
% This program accompanies the book:
%
%          C. Pozrikidis
% ``Numerical Computation in Science and Engineering''
%         Oxford University Press
%------------------------------------------------
 
%===============================================
% anination of a 2D line in the xy plane
% interpolated by a clamped spline
%
% LEGEND:
% ------
%
% N:        Number of intervals
% xp(i):    x value of ith point
% fp(i):    f value of ith point
% a, b, c:  cubic spline coefficients
%
% kpl:  number of intervals
%       for plotting the spline
%===============================================

%=====
% point coordinates
%=====

% N+1 points:

fpx=[0.00 0.12  0.27  0.32  0.45  0.55  0.68  0.71  0.80  0.87 1.01 1.5];
fpy=[0.07 0.22  0.43  0.34  0.49  0.54  0.49  0.32  0.15  0.11 0.10 0.99];
fpy=[0.00 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0 0.0 0.0];

N = 11;

%---
% initialize the point velocity
%---

for i=1:N+1
 velx(i) = 1.0;
 vely(i) = 0.9;
end

%---
% animation loop
%---

for repeat=1:9000

%---
% points move randomly
%---

Ido = 0;

if(Ido==1)
 Dt = 0.01;
 for i=1:N+1
 fpx(i) = fpx(i)+Dt*(rand(1)-0.5);
 fpy(i) = fpy(i)+Dt*(rand(1)-0.5);
 end
end

%---
% points bounce inside a box
%---

Ido =0;

 if(Ido==1)
  Dt = 0.01;
  for i=1:N+1
   fpx(i) = fpx(i)+velx(i)*Dt;
   fpy(i) = fpy(i)+vely(i)*Dt;
   if(fpx(i)> win) velx(i)=-abs(velx(i)); end
   if(fpx(i)<  0) velx(i)= abs(velx(i)); end
   if(fpy(i)> win) vely(i)=-abs(vely(i)); end
   if(fpy(i)<  0) vely(i)= abs(vely(i)); end
  end
end

%---
% point source
%---

Ido = 1;

 if(Ido==1)
  Dt = 0.0001;
  xsrc = 0.0;
  ysrc = 1.0;
  for i=1:N+1
   R2 = (fpx(i)-xsrc)^2+(fpy(i)-ysrc)^2;
   velx(i) = -(fpx(i)-0.0)/R2;
   vely(i) = -(fpy(i)-1.0)/R2;
   fpx(i) = fpx(i)+velx(i)*Dt;
   fpy(i) = fpy(i)+vely(i)*Dt;
  end
end

win =1.5;
axis([0 win 0 ysrc])

%====
% compute the polygonal arc length
% to be used as the interpolation variable
%====
 
xp(1) = 0.0D0;
 
for i=2:N+1
 xp(i) = xp(i-1)+ sqrt((fpx(i)-fpx(i-1))^2+(fpy(i)-fpy(i-1))^2);
end

%====
% spline coefficients
%====

slope1 = 1.0; slope2=1.0;
[ax,bx,cx] = splc_clm (N,xp,fpx,slope1,slope2);
slope1 = 0.0; slope2=0.0;
[ay,by,cy] = splc_clm (N,xp,fpy,slope1,slope2);

%====
% plotting
%====

kpl=128; % number of plotting intervals

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

 xint = xint+dxint;

end

%---
% animation
%---

if(repeat==1) 
 Handle1 = plot(xx,yy,'-','linewidth',3);
 box
 hold on
 Handle2 = plot(fpx,fpy,'ro');
 set(Handle1,'EraseMode','xor')
 set(Handle2,'EraseMode','xor')
 xlabel('x')
 ylabel('y')
end

set(Handle1,'XData',xx,'YData',yy);
set(Handle2,'XData',fpx,'YData',fpy)
drawnow 

%--
end % of animation loop
%--
