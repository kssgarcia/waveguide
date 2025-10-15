function [Gxx,Gxy,Gxz ...
         ,Gyx,Gyy,Gyz ...
         ,Gzx,Gzy,Gzz ...
         ,px,py,pz ...
         ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz ...
         ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz ...
         ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz ...
         ] = ...
...
      sgf_3d_fs ...
...
       (Iopt ...
       ,x,y,z ...
       ,x0,y0,z0 ...
       )

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%---------------------------------------
% Free-space Green's function: Stokeslet
%
% Pozrikidis (1992, p. 23)
%
% Iopt =  1 generates only the Green's function
%      ne 1 generates the Green's function,
%           pressure, and stress
%---------------------------------------

%---
% initialize
%---

     px=0; py=0; pz=0;
     Txxx=0; Txxy=0; Txxz=0; Tyxy=0; Tyxz=0; Tzxz=0;
     Txyx=0; Txyy=0; Txyz=0; Tyyy=0; Tyyz=0; Tzyz=0;
     Txzx=0; Txzy=0; Txzz=0; Tyzy=0; Tyzz=0; Tzzz=0;

      dx = x-x0;
      dy = y-y0;
      dz = z-z0;

      dxx = dx*dx;
      dxy = dx*dy;
      dxz = dx*dz;
      dyy = dy*dy;
      dyz = dy*dz;
      dzz = dz*dz;

      r  = sqrt(dxx+dyy+dzz);
      r3 = r*r*r;

      ri  = 1.0D0/r;
      ri3 = 1.0D0/r3;

      Gxx = ri + dxx*ri3;
      Gxy =      dxy*ri3;
      Gxz =      dxz*ri3;
      Gyy = ri + dyy*ri3;
      Gyz =      dyz*ri3;
      Gzz = ri + dzz*ri3;

      Gyx = Gxy;
      Gzx = Gxz;
      Gzy = Gyz;

%--------------
% stress tensor
%--------------

   if(Iopt==2) 

      cf = -6.0D0/r^5;

      Txxx = dxx*dx * cf;
      Txxy = dxy*dx * cf;
      Txxz = dxz*dx * cf;
      Tyxy = dyy*dx * cf;
      Tyxz = dyz*dx * cf;
      Tzxz = dzz*dx * cf;

      Txyx = Txxy;
      Txyy = Tyxy;
      Txyz = Tyxz;
      Tyyy = dyy*dy * cf;
      Tyyz = dyz*dy * cf;
      Tzyz = dzz*dy * cf;

      Txzx = Txxz;
      Txzy = Tyxz;
      Txzz = Tzxz;
      Tyzy = dyy*dz * cf;
      Tyzz = dyz*dz * cf;
      Tzzz = dzz*dz * cf;

%---------
% pressure
%---------

      cf = 2.0D0*ri3;
      px = dx * cf;
      py = dy * cf;
      pz = dz * cf;

   end

%-----
% done
%-----

   return
