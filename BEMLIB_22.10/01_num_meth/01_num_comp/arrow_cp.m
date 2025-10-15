function gaidaros = arrow_cp (x1,y1,dx,dy)

%========================================
% FDLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved
%
% This program is to be used only under the
% stipulations of the licensing agreement
%========================================

%=====================================================
% Generate coordinates for plotting a five-point arrow
% in the plane starting at the point(x1,x2)
% and ending at the point (x1+dx, y1+dy)
%
% SYMBOLS:
% --------
%
% dx, dy: arrow vector
%
% gaidaros(i,1): x-coordinate of the ith point
% gaidaros(i,2): y-coordinate of the ith point
%              where i=1,2,3,4,5
%
% angle:  angle of arrow tip in radians
%
% tip: length of arrow tip sides
%      as a fraction of the arrow length
%=====================================================

  angle = 0.3;
  tip = 0.50;

%--------
% prepare
%--------

  x2 = x1+dx;
  y2 = y1+dy;
  cs = cos(angle);
  sn = sin(angle);
  dxi = -dx;
  dyi = -dy;

%------------
% first point
%------------

  gaidaros(1,1) = x1;
  gaidaros(1,2) = y1;

%---
% second point
%---

  gaidaros(2,1) = x2;
  gaidaros(2,2) = y2;

%---
% third point
%---

  gaidaros(3,1) = x2+( dxi*cs+dyi*sn)*tip;
  gaidaros(3,2) = y2+(-dxi*sn+dyi*cs)*tip;

%---
% fourth point
%---

  gaidaros(4,1) = x2;
  gaidaros(4,2) = y2;

%---
% fifth point
%---

  gaidaros(5,1) = x2+(dxi*cs-dyi*sn)*tip;
  gaidaros(5,2) = y2+(dxi*sn+dyi*cs)*tip;

%-----
% done
%-----

  return
