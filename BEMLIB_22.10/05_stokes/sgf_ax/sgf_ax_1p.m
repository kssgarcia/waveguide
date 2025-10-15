function [GXX,GXY,GYX,GYY ...
...
  ,QXXX,QXXY,QXYX,QXYY ...
  ,QYXX,QYXY,QYYX,QYYY ...
  ...
  ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
...
  = sgf_ax_1p (Iopt ...
              ,X0,Y0 ...
              ,X,Y ...
              ,L ...
              ,Nsum,Np)

%-----------------------------------------
% FDLIB, BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement.
%----------------------------------------

%======================================================
%  Axisymmetric periodic Green's function of Stokes flow
%
%  Symbols:
%  -------
%
%  L: 		period along the x axis
%  Nsum, Np: 	summation range limits
%
%  y stands for radial distance "sigma"
%
%-------------
%
%  Let (bx, by) be the strength of the point-force ring.
%  The induced velocity is:
%
%  ux(x0) = SXX(x0,x) * bx + SXY(x0,x) * by
%  uy(x0) = SYX(x0,x) * bx + SYY(x0,x) * by
%
%  The kernel of the axisymmetric double-layer potential is:
%
%   Idlpx(x0) = ux * ( Qxxx * vnx + Qxxy * vny)
%             + uy * ( Qxyx * vnx + Qxyy * vny)
%
%   Idlpy(x0) = ux * ( Qyxx * vnx + Qyxy * vny)
%             + uy * ( Qyyx * vnx + Qyyy * vny)
%
%   This is the flow due to a ring distribution of stresslets
%
%   Summation will be done over x
%
%------------
%
%  Iopt = 1 only the Green's function
%  Iopt = 2          Green's function and the stress tensor
%======================================================

%----------
% constants
%----------

  pi4 = 4.0*pi;

%--------
% prepare
%--------

 subtr = 0.0;
 subtr = pi4*Y;

%---------------------------
% prepare for fast summation
% using Aitken extrapolation
%---------------------------

 N0 = Nsum;
 N1 = N0*Np;
 N2 = N1*Np;

%-------------
% primary ring
%-------------

  [GXX,GXY ...
  ,GYX,GYY ...
  ,QXXX,QXXY,QXYX,QXYY ...
  ,QYXX,QYXY,QYYX,QYYY ...
  ,PXX,PXY,PYX,PYY ...
  ,Iaxis] ...
  ...
   = sgf_ax_fs (Iopt ...
               ,X0,Y0 ...
               ,X,Y);

%   GXX = GXX - subtr/abs(X+L);

%------------------------
% sum over periodic rings
%------------------------

  Xorg = X;

  for Ir=1:N2

   X = Xorg+Ir*L;     % ring on the right

   [SXX,SXY,SYX,SYY ...
   ,RXXX,RXXY,RXYX,RXYY ...
   ,RYXX,RYXY,RYYX,RYYY ...
   ,EXX,EXY,EYX,EYY ...
   ,Iaxis] ...
 ...
   = sgf_ax_fs (Iopt ...
               ,X0,Y0 ...
               ,X,Y);

%-------
       if(Ir==N0 |Ir==N1 | Ir==N2)  % hold for extrapolation

       FXX = SXX;
       FXY = SXY;
       FYX = SYX;
       FYY = SYY;

        if(Iopt==2)   % stress and pressure
         FXXX = RXXX;
         FXXY = RXXY;
         FXYX = RXYX;
         FXYY = RXYY;
         FYXX = RYXX;
         FYXY = RYXY;
         FYYX = RYYX;
         FYYY = RYYY;
         UXX  = EXX;
         UXY  = EXY;
         UYX  = EYX;
         UYY  = EYY;
        end

      end
%-------

      GXX = GXX + SXX;
      GXY = GXY + SXY;
      GYX = GYX + SYX;
      GYY = GYY + SYY;

      if(Iopt==2)

       QXXX = QXXX + RXXX;
       QXXY = QXXY + RXXY;
       QXYX = QXYX + RXYX;
       QXYY = QXYY + RXYY;
       QYXX = QYXX + RYXX;
       QYXY = QYXY + RYXY;
       QYYX = QYYX + RYYX;
       QYYY = QYYY + RYYY;

       PXX  = PXX+EXX;
       PXY  = PXY+EXY;
       PYX  = PYX+EYX;
       PYY  = PYY+EYY;

      end

      X = Xorg-Ir*L;    % ring on the left

     [SXX,SXY,SYX,SYY ...
     ,RXXX,RXXY,RXYX,RXYY ...
     ,RYXX,RYXY,RYYX,RYYY ...
     ,EXX,EXY,EYX,EYY ...
     ,Iaxis] ...
    ...
     = sgf_ax_fs (Iopt ...
                 ,X0,Y0 ...
                 ,X,Y);

%---
%   if(Ir==N2)
%    SXX
%   subtr/abs(X)
%   end
%   SXX = SXX - subtr/abs(X);
%---

% factor of 2 from left and righ ring:

    SXX = SXX - 2.0*subtr/abs(Ir*L);


%============
      if(Ir==N0 |Ir==N1 | Ir==N2)  % hold for extrapolation

       FXX = FXX + SXX;
       FXY = FXY + SXY;
       FYX = FYX + SYX;
       FYY = FYY + SYY;

       if(Iopt==2)
        FXXX = FXXX + RXXX;
        FXXY = FXXY + RXXY;
        FXYX = FXYX + RXYX;
        FXYY = FXYY + RXYY;
        FYXX = FYXX + RYXX;
        FYXY = FYXY + RYXY;
        FYYX = FYYX + RYYX;
        FYYY = FYYY + RYYY;
        UXX  = UXX  + EXX;
        UXY  = UXY  + EXY;
        UYX  = UYX  + EYX;
        UYY  = UYY  + EYY;
       end

      end
%============

      GXX = GXX + SXX;
      GXY = GXY + SXY;
      GYX = GYX + SYX;
      GYY = GYY + SYY;

      if(Iopt==2)
       QXXX = QXXX + RXXX;
       QXXY = QXXY + RXXY;
       QXYX = QXYX + RXYX;
       QXYY = QXYY + RXYY;
       QYXX = QYXX + RYXX;
       QYXY = QYXY + RYXY;
       QYYX = QYYX + RYYX;
       QYYY = QYYY + RYYY;
       PXX  = PXX+EXX;
       PXY  = PXY+EXY;
       PYX  = PYX+EYX;
       PYY  = PYY+EYY;
      end

      if(Ir==N0)
        GXX0 = GXX - 0.5D0 * FXX;
        GXY0 = GXY - 0.5D0 * FXY;
        GYX0 = GYX - 0.5D0 * FYX;
        GYY0 = GYY - 0.5D0 * FYY;
      elseif(Ir==N1)
        GXX1 = GXX - 0.5D0 * FXX;
        GXY1 = GXY - 0.5D0 * FXY;
        GYX1 = GYX - 0.5D0 * FYX;
        GYY1 = GYY - 0.5D0 * FYY;
      elseif(Ir==N2)
        GXX2 = GXX - 0.5D0 * FXX;
        GXY2 = GXY - 0.5D0 * FXY;
        GYX2 = GYX - 0.5D0 * FYX;
        GYY2 = GYY - 0.5D0 * FYY;
      end

%---
      if(Iopt==2)
       if(Ir==N0)
        QXXX0 = QXXX - 0.5D0 * FXXX;
        QXXY0 = QXXY - 0.5D0 * FXXY;
        QXYX0 = QXYX - 0.5D0 * FXYX;
        QXYY0 = QXYY - 0.5D0 * FXYY;
        QYXX0 = QYXX - 0.5D0 * FYXX;
        QYXY0 = QYXY - 0.5D0 * FYXY;
        QYYX0 = QYYX - 0.5D0 * FYYX;
        QYYY0 = QYYY - 0.5D0 * FYYY;
         PXX0 = PXX  - 0.5D0 * UXX;
         PXY0 = PXY  - 0.5D0 * UXY;
         PYX0 = PYX  - 0.5D0 * UYX;
         PYY0 = PYY  - 0.5D0 * UYY;
       elseif(Ir==N1)
        QXXX1 = QXXX - 0.5D0 * FXXX;
        QXXY1 = QXXY - 0.5D0 * FXXY;
        QXYX1 = QXYX - 0.5D0 * FXYX;
        QXYY1 = QXYY - 0.5D0 * FXYY;
        QYXX1 = QYXX - 0.5D0 * FYXX;
        QYXY1 = QYXY - 0.5D0 * FYXY;
        QYYX1 = QYYX - 0.5D0 * FYYX;
        QYYY1 = QYYY - 0.5D0 * FYYY;
         PXX1 = PXX  - 0.5D0 * UXX;
         PXY1 = PXY  - 0.5D0 * UXY;
         PYX1 = PYX  - 0.5D0 * UYX;
         PYY1 = PYY  - 0.5D0 * UYY;
       elseif(Ir==N2)
        QXXX2 = QXXX - 0.5D0 *FXXX;
        QXXY2 = QXXY - 0.5D0 *FXXY;
        QXYX2 = QXYX - 0.5D0 *FXYX;
        QXYY2 = QXYY - 0.5D0 *FXYY;
        QYXX2 = QYXX - 0.5D0 *FYXX;
        QYXY2 = QYXY - 0.5D0 *FYXY;
        QYYX2 = QYYX - 0.5D0 *FYYX;
        QYYY2 = QYYY - 0.5D0 *FYYY;
         PXX2 = PXX  - 0.5D0 *UXX;
         PXY2 = PXY  - 0.5D0 *UXY;
         PYX2 = PYX  - 0.5D0 *UYX;
         PYY2 = PYY  - 0.5D0 *UYY;
       end
      end
%---

   end

%---
% Aitken extrapolation for the velocity
%---

 GXX = (GXX0*GXX2-GXX1^2)/(GXX2-2.0D0*GXX1+GXX0);
 GXY = (GXY0*GXY2-GXY1^2)/(GXY2-2.0D0*GXY1+GXY0);
 GYX = (GYX0*GYX2-GYX1^2)/(GYX2-2.0D0*GYX1+GYX0);
 GYY = (GYY0*GYY2-GYY1^2)/(GYY2-2.0D0*GYY1+GYY0);

%---
% Aitken extrapolation for the stress and pressure
%---

 if(Iopt==2)

  QXXX = (QXXX0*QXXX2-QXXX1^2)/(QXXX2-2.0D0 *QXXX1+QXXX0);
  QXXY = (QXXY0*QXXY2-QXXY1^2)/(QXXY2-2.0D0 *QXXY1+QXXY0);
  QXYX = (QXYX0*QXYX2-QXYX1^2)/(QXYX2-2.0D0 *QXYX1+QXYX0);
  QXYY = (QXYY0*QXYY2-QXYY1^2)/(QXYY2-2.0D0 *QXYY1+QXYY0);

  QYXX = (QYXX0*QYXX2-QYXX1^2)/(QYXX2-2.0D0 *QYXX1+QYXX0);
  QYXY = (QYXY0*QYXY2-QYXY1^2)/(QYXY2-2.0D0 *QYXY1+QYXY0);
  QYYX = (QYYX0*QYYX2-QYYX1^2)/(QYYX2-2.0D0 *QYYX1+QYYX0);
  QYYY = (QYYY0*QYYY2-QYYY1^2)/(QYYY2-2.0D0 *QYYY1+QYYY0);

  PXX = (PXX0*PXX2-PXX1^2)/(PXX2-2.0D0 *PXX1+PXX0);
  PXY = (PXY0*PXY2-PXY1^2)/(PXY2-2.0D0 *PXY1+PXY0);
  PYX = (PYX0*PYX2-PYX1^2)/(PYX2-2.0D0 *PYX1+PYX0);
  PYY = (PYY0*PYY2-PYY1^2)/(PYY2-2.0D0 *PYY1+PYY0);

 end

%---
% restore
%---

   X = Xorg;

%-----
% done
%-----

return
