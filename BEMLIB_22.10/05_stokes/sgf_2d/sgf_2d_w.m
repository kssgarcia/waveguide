function [Gxx,Gxy ...
         ,Gyx,Gyy ...
         ,Px,Py ...
         ,Txxx,Txxy,Tyxx,Tyxy ...
         ,Txyx,Txyy,Tyyx,Tyyy ] ...
 ...
  = sgf_2d_w  ...
       ...
      (x,y ...
      ,x0,y0 ...
      ,wall ...
      ,Iopt ...
      )

%-----------------------------------------
% FDLIB BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%--------------------------------------------------
% Green's function for flow bounded by a plane wall
% located at y = wall
%
% see Pozrikidis (1992, p. 93)
%
% Iopt =  1 compute only G
%      ne 1 compute G, P, and T
%--------------------------------------------------

%---
% initialize
%---

Px=0; Py=0;

Txxx=0;
Txxy=0;
Tyxx=0;
Tyxy=0;
Txyx=0;
Txyy=0;
Tyyx=0;
Tyyy=0;

%====================
% primary point force
%====================

 dx = x-x0;
 dy = y-y0;
 dxx = dx*dx;
 dxy = dx*dy;
 dyy = dy*dy;
 r2  = dxx+dyy;
 r   = sqrt(r2);
 ri2 = 1.0D0/r2;
 rl = log(r);

 Gxx = -rl + dxx*ri2;
 Gxy =       dxy*ri2;
 Gyy = -rl + dyy*ri2;

%--------------------
% pressure and stress:
%--------------------

if(Iopt~=1)
 cf = -4.0D0*ri2*ri2;
 Txxx = dxx*dx*cf;
 Txxy = dxy*dx*cf;
 Tyxy = dyy*dx*cf;
 Txyx = Txxy;
 Txyy = Tyxy;
 Tyyy = dyy*dy*cf;
 Px = 2.0D0*dx*cf*ri2;
 Py = 2.0D0*dy*cf*ri2;
end

%==================
% image point force
%==================

 y0im = 2.0D0*wall-y0;
 dy = y-y0im;
 dxy = dx*dy;
 dyy = dy*dy;
 r2 = dxx+dyy;
 r  = sqrt(r2);
 ri2 = 1.0D0/r2;
 ri4 = ri2*ri2;
 rl = log(r);

 Gxx = Gxx + rl - dxx*ri2;
 Gxy = Gxy      - dxy*ri2;
 Gyx = Gxy;
 Gyy = Gyy + rl - dyy*ri2;

%--------------------
% pressure and stress:
%--------------------

if(Iopt~=1)

cf = -4.0D0*ri4;

Txxx = Txxx - dxx*dx*cf;
Txxy = Txxy - dxy*dx*cf;
Tyxy = Tyxy - dyy*dx*cf;
Txyx = Txxy;
Txyy = Tyxy;
Tyyy = Tyyy - dyy*dy*cf;

Px = Px-2.0D0*dx*cf*ri2;
Py = Py-2.0D0*dy*cf*ri2;

end

%=======================
% image potential dipole
%=======================

cf = 2.0D0*ri4;

DPxx =   ri2 - dxx*cf;
DPyx =       - dxy*cf;

DPxy = - DPyx;
DPyy = - ri2 + dyy*cf;

%--------------------
% pressure and stress:
%--------------------

if(Iopt~=1)

 cf  =  4.0D0*ri4;
 cf1 = 16.0D0*ri2*ri2*ri2;

 TDPxxx = - 3.0D0*dx*cf + dxx*dx*cf1;
 TDPxxy = - dy      *cf + dxy*dx*cf1;
 TDPyxy = - dx      *cf + dyy*dx*cf1;

 TDPxyx =   dy      *cf - dxx*dy*cf1;
 TDPxyy =   dx      *cf - dxy*dy*cf1;
 TDPyyy =   3.0D0*dy*cf - dyy*dy*cf1;

end

%=======================
% image Stokeslet dipole
%=======================

SDxx = dy*DPxx;
SDyx = dy*DPyx-dx*ri2;

SDxy = dy*DPxy-dx*ri2;
SDyy = dy*DPyy;

%--------------------
% pressure and stress:
%--------------------

if(Iopt~=1)

  cf = 4.0D0*ri4;

  TSDxxx = dy * TDPxxx + dxy * cf;
  TSDxxy = dy * TDPxxy;
  TSDyxy = dy * TDPyxy + dxy * cf;

  TSDxyx = dy * TDPxyx - (dyy-dxx)*cf;
  TSDxyy = dy * TDPxyy - (   -dxy)*cf;
  TSDyyy = dy * TDPyyy;

  PSDx =      - 4.0D0*dxy*ri4;
  PSDy = -2.0D0*(dxx-dyy)*ri4;

end

%====================
% putting it together
%====================

h0   = y0-wall;
h02  = 2.0D0*h0;
h0s2 = 2.0D0*h0*h0;

Gxx = Gxx + h0s2*DPxx - h02*SDxx;
Gxy = Gxy + h0s2*DPxy - h02*SDxy;

Gyx = Gyx + h0s2*DPyx - h02*SDyx;
Gyy = Gyy + h0s2*DPyy - h02*SDyy;

%--------------------
% pressure and stress
%--------------------

if(Iopt~=1)

Txxx = Txxx  ...
     + h0s2 * TDPxxx  ...
     - h02  * TSDxxx;
Txxy = Txxy  ...
     + h0s2 * TDPxxy  ...
     - h02  * TSDxxy;
Tyxy = Tyxy  ...
     + h0s2 * TDPyxy  ...
     - h02  * TSDyxy;
Txyx = Txyx  ...
     + h0s2 * TDPxyx  ...
     - h02  * TSDxyx;
Txyy = Txyy  ...
     + h0s2 * TDPxyy  ...
     - h02  * TSDxyy;
Tyyy = Tyyy  ...
     + h0s2 * TDPyyy  ...
     - h02  * TSDyyy;
Tyxx = Txxy;
Tyyx = Txyy;

Px = Px - h02*PSDx;
Py = Py - h02*PSDy;

end

%-----
% done
%-----

return
