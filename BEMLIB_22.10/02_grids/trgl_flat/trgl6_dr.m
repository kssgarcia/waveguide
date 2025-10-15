clear all
close all

%=========================================
% FDLIB  BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%=========================================

%=========================================
% Driver for triangulating:
%
% 1) a square
% 2) a square with a square hole
% 3) a square with a circular hole
% 4) a disk
%=========================================

ishape = 2;  % square with a square hole
ishape = 1;  % square
ishape = 3;  % square with a circular hole
ishape = 4; % disk

%---
% graphics
%---

figure(1)
hold on

%---
% parameters
%---

a = 0.3;
ndiv=1;

%---
% triangulation
%---

if(ishape==1)
 [ne,ng,p,c,efl,gfl] = trgl6_sqr(ndiv);
elseif(ishape==2)
 [ne,ng,p,c,efl,gfl] = trgl6_ss(a,ndiv)
elseif(ishape==3)
 [ne,ng,p,c,efl,gfl] = trgl6_sc(a,ndiv);
elseif(ishape==4)
 [ne,ng,p,c,efl,gfl] = trgl6_disk(ndiv);
end

%---
% warp into 3D
%---

for i=1:ng
 p(i,3) = 0.0;
 if(ishape==2)
  rrr = sqrt(p(i,1)^2+p(i,2)^2);
  p(i,3) = 0.5*a/rrr^2;
  if(gfl(i)==1) p(i,3) = 0.0; end
  if(gfl(i)==2) p(i,3) = 0.5/a; end
 end
end

%---
% run over elements
% compute alpha, beta, gamma
%---

  for i=1:ne
   i1 = c(i,1);
   i2 = c(i,2);
   i3 = c(i,3);
   i4 = c(i,4);
   i5 = c(i,5);
   i6 = c(i,6);
   xp(1)=p(i1,1);
   xp(2)=p(i4,1);
   xp(3)=p(i2,1);
   xp(4)=p(i5,1);
   xp(5)=p(i3,1);
   xp(6)=p(i6,1);
   xp(7)=p(i1,1);
   yp(1)=p(i1,2);
   yp(2)=p(i4,2);
   yp(3)=p(i2,2);
   yp(4)=p(i5,2);
   yp(5)=p(i3,2);
   yp(6)=p(i6,2);
   yp(7)=p(i1,2);
   zp(1)=p(i1,3);
   zp(2)=p(i4,3);
   zp(3)=p(i2,3);
   zp(4)=p(i5,3);
   zp(5)=p(i3,3);
   zp(6)=p(i6,3);
   zp(7)=p(i1,3);

   [al(i), be(i), ga(i)] = elm6_2d_abc ...
...
    (xp(1),yp(1),zp(1), xp(2),yp(2),zp(2), xp(3),yp(3),zp(3) ...
    ,xp(4),yp(4),zp(4) ...
    ,xp(5),yp(5),zp(5), xp(6),yp(6),zp(6));

   plot3(xp,yp,zp,'k');
   hold on;
   plot3(xp,yp,zp,'ko');

  end

  xlabel('x');
  ylabel('y');
  axis equal
