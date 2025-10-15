function [Gxx,Gxy ...
         ,Gyx,Gyy ...
         ,Px,Py ...
         ,TXXX,TXXY,TYXX,TYXY ...
         ,TXYX,TXYY,TYYX,TYYY ] ...
 ...
   = sgf_2d_1p_w  ...
   ...
   (X,Y ...
   ,X0,Y0 ...
   ,wall ...
   ,period...
   ,Iopt ...
   )

%-----------------------------------------
% FDLIB BEMLIB
%
% Copyright by C. Pozrikidis, 1999
% All rights reserved.
%
% This program is to be used only under the
% stipulations of the licensing agreement
%----------------------------------------

%------------------------------------------------- 
%  Green's function of two-dimensional Stokes flow
%  associated with an array of point forces 
%  deployed along the x axis above a plane wall
%  located at y = wall
%
%  One of the point forces is located at (X0,Y0)
%
%  The field point is located at (X,Y)
%
%  SYMBOLS:
%  -------
%
%  period:  Distance between the point forces
%
%  Iselect = 1 computes only G
%           ne computes G, p, T
%
%  Ising = 1 computes G-S where S is the Stokeslet
%            similarly for the stress and pressure 
%------------------------------------------------- 

%---
% initialize
%---

Px = 0.0; Py = 0.0;

TXXX = 0; TXXY = 0;
TYXX = 0; TYXY = 0;

TXYX = 0; TXYY = 0;
TYYX = 0; TYYY = 0;

%------------
% wave number
%------------

 wn = 2*pi/period;

%-------------------------------
% periodic array of point forces
%-------------------------------

 DX = X-X0;
 DY = Y-Y0;

 CHY = cosh(wn*DY);
 SHY = sinh(wn*DY);
 CX  =  cos(wn*DX);
 SX  =  sin(wn*DX);

 D   = CHY-CX;

 A  = 0.5D0*log(2.0D0*D);

 AX = 0.5D0*wn*SX /D;
 AY = 0.5D0*wn*SHY/D;

 Gxx = -A-DY*AY+1.0D0;
 Gxy =  DY*AX;
 Gyx =  Gxy;
 Gyy = -A+DY*AY;

%------------------------------------
% compute the stress and the pressure
%------------------------------------

 if(Iopt>1)

  D2  = D*D;
  AYY = 0.5D0*wn*wn*(1.0D0-CX*CHY)/D2;
  AXY =-0.5D0*wn*wn*       SX*SHY /D2;

  T1 = -2.0D0*(2.0D0*AX+DY*AXY);
  T2 = -2.0D0*(AY+DY*AYY);
  T3 =  2.0D0*DY*AXY;
  T4 = -2.0D0*(AY-DY*AYY);

  Px = 2.0D0*AX;
  Py = 2.0D0*AY;
 
end

%--------------------
% image singularities
%--------------------

 YW = Y+Y0-2.0D0*wall;

 CHW = cosh(wn*YW);
 SHW = sinh(wn*YW);
 DW  = CHW-CX;
 DW2 = DW*DW;
 DW3 = DW2*DW;

 AW   =  0.5D0*log(2.0D0*DW);

 AWX  =  0.5D0*wn*SX /DW;
 AWY  =  0.5D0*wn*SHW/DW;
 AWYY =  0.5D0*wn*wn*(1.0D0-CX*CHW)/DW2;
 AWXY = -0.5D0*wn*wn*       SX*SHW /DW2;

 h0  = Y0-wall;
 h02 = h0*h0;

 Gxx = Gxx + AW-1.0D0+YW*AWY ...
       +2.0D0*h0*YW*AWYY ...
       -2.0D0*h02*AWYY;
 Gxy = Gxy - YW*AWX ...
       +2.0D0*h0*(AWX+YW*AWXY) ...
       -2.0D0*h02*AWXY;
 Gyx = Gyx - YW*AWX ...
       +2.0D0*h0*(AWX-YW*AWXY) ...
       +2.0D0*h02*AWXY;
 Gyy = Gyy + AW-YW*AWY ...
       +2.0D0*h0*YW*AWYY ...
       -2.0D0*h02*AWYY;

%---
 if(Iopt>1)
%---

 cf = wn^3/(2.0D0*DW3);

 AWXYY= SX *(CHW^2+CHW*CX-2.0D0)*cf;
 AWXXY=-SHW*( CX^2+CHW*CX-2.0D0)*cf;

 PWX = 2.0D0*(-AWX-2.0D0*h0*AWXY);
 PWY = 2.0D0*(-AWY+2.0D0*h0*AWYY);

     SXXX = AWX+YW*AWXY ...
          +2.0D0*h0*YW*AWXYY ...
          -2.0D0*h02*AWXYY;
     SXYX = YW*AWYY ...
          -2.0D0*h0*(AWYY-YW*AWXXY) ...
          -2.0D0*h02*AWXXY;
     SYXX = YW*AWYY ...
          -2.0D0*h0*(AWYY+YW*AWXXY) ...
          +2.0D0*h02*AWXXY;
     SYYX = AWX-YW*AWXY ...
          +2.0D0*h0*YW*AWXYY ...
          -2.0D0*h02*AWXYY;
     SXXY = 2.0D0*AWY+YW*AWYY ...
          +2.0D0*h0*(AWYY-YW*AWXXY) ...
          +2.0D0*h02*AWXXY;
     SXYY = -AWX-YW*AWXY ...
          +2.0D0*h0*(2.0D0*AWXY+YW*AWXYY) ...
          -2.0D0*h02*AWXYY;
     SYXY = -AWX-YW*AWXY ...
          -2.0D0*h0*YW*AWXYY ...
          +2.0D0*h02*AWXYY;
     SYYY=       -YW*AWYY ...
          +2.0D0*h0*(AWYY-YW*AWXXY) ...
          +2.0D0*h02*AWXXY;

      TXXX = - PWX + 2.0D0*SXXX;
      TXXY =        SXXY+SYXX;
      TYXX =   TXXY;
      TYXY = - PWX + 2.0D0*SYXY;

      TXYX = - PWY + 2.0D0*SXYX;
      TXYY =        SXYY+SYYX;
      TYYX =   TXYY;
      TYYY = - PWY + 2.0D0*SYYY;

      TXXX = TXXX + T1;
      TXXY = TXXY + T2;
      TYXX = TYXX + T2;
      TYXY = TYXY + T3;

      TXYX = TXYX + T2;
      TXYY = TXYY + T3;
      TYYX = TYYX + T3;
      TYYY = TYYY + T4;

      Px = Px+PWX;
      Py = Py+PWY;

%---
end
%---

%-----
% done
%-----

return
