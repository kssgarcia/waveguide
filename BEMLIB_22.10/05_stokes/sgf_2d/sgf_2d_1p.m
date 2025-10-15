function [Gxx,Gxy ...
         ,Gyx,Gyy ...
         ,Px,Py ...
         ,Txxx,Txxy,Tyxx,Tyxy ...
         ,Txyx,Txyy,Tyyx,Tyyy ] ...
 ...
  = sgf_2d_1p  ...
    ...
    (X,Y ...
    ,X0,Y0 ...
    ,period ...
    ,Iopt ...
    ,Idesing ...
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

%------------------------------------------------- 
%  Green's function of two-dimensional Stokes flow
%  associated with an array
%  of point forces arranged along the x axis
%
%  One point force is located at: (X0,Y0)
%
%  The field point is located at: (X,Y)
%
%  SYMBOLS:
%  -------
%
%  period: distance between two point forces
%
%  Iopt = 1  Computes only G
%       ne 1  Computes G, p, T
%
%  If Idesing = 1 will compute G-S,
%                 where S is the Stokeslet
%                 similarly for the stress and pressure
%------------------------------------------------- 

%------------
% wave number
%------------

 wn = 2.0*pi/period;

%----------
% compute G
%----------

  DX = X-X0;
  DY = Y-Y0;

  CHY = cosh(wn*DY);
  SHY = sinh(wn*DY);
  CX  =  cos(wn*DX);
  SX  =  sin(wn*DX);

  D = CHY-CX;

  A  = 0.5D0*log(2.0D0*D);
  AX = 0.5D0*wn*SX /D;
  AY = 0.5D0*wn*SHY/D;

  Gxx = -A-DY*AY;   % +1.0D0  % optional
  Gxy =  DY*AX;
  Gyx =  Gxy;
  Gyy = -A+DY*AY;

%---
% subtract out the Stokeslet
% if desired
%---

 if(Idesing==1) 

   DX2   = DX*DX;
   DY2   = DY*DY;
   DXY   = DX*DY;
   DIST2 = DX2 + DY2;
   SING  = 0.5D0*log(DIST2);

   Gxx = Gxx + SING - DX2/DIST2;
   Gxy = Gxy        - DXY/DIST2;
   Gyx = Gyx        - DXY/DIST2;
   Gyy = Gyy + SING - DY2/DIST2;

 end

%---------------------------------
%  compute the stress and pressure
%  tensors
%---------------------------------

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

%---
%  subtract out the Stokeslet
%---

      if(Idesing==1)

          cf  = 4.0D0/DIST2^2;
          T1 = T1 + cf*DX*DX*DX;
          T2 = T2 + cf*DX*DX*DY;
          T3 = T3 + cf*DX*DY*DY;
          T4 = T4 + cf*DY*DY*DY;

          Px = Px - 2.0D0*DX/DIST2;
          Py = Py - 2.0D0*DY/DIST2;

      end

%-------
% finish
%-------

      Txxx = T1;
      Txxy = T2;
      Tyxx = T2;
      Tyxy = T3;

      Txyx = T2;
      Txyy = T3;
      Tyyx = T3;
      Tyyy = T4;

      end

%-----
% done
%-----

   return
