      subroutine sgf_2d_1p_w 
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,wall
     +   ,RL
     +   ,Idesing
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,Px,Py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------------- 
c  Green's function of two-dimensional Stokes flow
c  associated with an array of point forces 
c  deployed along the x axis above a plane wall
c  located at y = wall
c
c  One of the point forces is located at: (X0,Y0)
c
c  The field point is located at: (X,Y)
c
c  SYMBOLS:
c  -------
c
c  RL:  Distance between the point forces
c
c  Iopt = 1 compute only G
c      ne 1 compute G, p, T
c
c  Idesing = 1 computes G-S where S is the Stokeslet
c              similarly for the stress and pressure 
c------------------------------------------------- 

      Implicit Double Precision (a-h,o-z)
 
c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

c------------
c wave number
c------------

      wn = pi2/RL

c-------------------------------
c periodic array of point forces
c-------------------------------

      DX = X-X0
      DY = Y-Y0

      CHY = cosh(wn*DY)
      SHY = sinh(wn*DY)
      CX  =  cos(wn*DX)
      SX  =  sin(wn*DX)

      D   = CHY-CX 

      A  = 0.5D0*Dlog(2.0D0*D)

      AX = 0.5D0*wn*SX /D
      AY = 0.5D0*wn*SHY/D

      Gxx = -A-DY*AY+1.0D0
      Gxy =  DY*AX
      Gyx =  Gxy
      Gyy = -A+DY*AY

c---
c subtract out the Stokeslet
c---

      if(Idesing.eq.1) then

        DX2   = DX*DX
        DY2   = DY*DY
        DXY   = DX*DY
        DIST2 = DX2+DY2
        SING  = 0.5D0*Dlog(DIST2)

        Gxx = Gxx + SING - DX2/DIST2
        Gxy = Gxy        - DXY/DIST2
        Gyx = Gyx        - DXY/DIST2 
        Gyy = Gyy + SING - DY2/DIST2

      end if

c----------------------------
c compute stress and pressure
c----------------------------

      if(Iopt.gt.1) then

       D2  = D*D
       AYY = 0.5D0*wn*wn*(1.0D0-CX*CHY)/D2
       AXY =-0.5D0*wn*wn*       SX*SHY /D2

       T1 = -2.0D0*(2.0D0*AX+DY*AXY)
       T2 = -2.0D0*(AY+DY*AYY)
       T3 =  2.0D0*DY*AXY
       T4 = -2.0D0*(AY-DY*AYY)

       Px = 2.0D0*AX
       Py = 2.0D0*AY

c----
c subtract out the Stokeslet
c----

      if(Idesing.eq.1) then

          cf = 4.0/DIST2**2
          T1 = T1 + cf*DX*DX*DX
          T2 = T2 + cf*DX*DX*DY
          T3 = T3 + cf*DX   *DY*DY
          T4 = T4 + cf      *DY*DY*DY

          Px = Px - 2.0D0*DX/DIST2
          Py = Py - 2.0D0*DY/DIST2

        end if

c-----------
      end if
c-----------

c--------------------
c image singularities
c--------------------

      YW = Y+Y0-2.0D0*wall

      CHW = Dcosh(wn*YW)
      SHW = Dsinh(wn*YW)
      DW  = CHW-CX
      DW2 = DW**2
      DW3 = DW2*DW

      AW   =  0.5D0*LOG(2.0D0*DW)

      AWX  =  0.5D0*wn*SX /DW
      AWY  =  0.5D0*wn*SHW/DW
      AWYY =  0.5D0*wn*wn*(1.0D0-CX*CHW)/DW2
      AWXY = -0.5D0*wn*wn*       SX*SHW /DW2

      h0  = Y0-wall
      h02 = h0*h0

      Gxx = Gxx + AW-1.0D0+YW*AWY
     +          -2.0D0*h02*AWYY
     +          +2.0D0*h0*YW*AWYY
      Gxy = Gxy - YW*AWX
     +          -2.0D0*h02*AWXY
     +          +2.0D0*h0*(AWX+YW*AWXY)
      Gyx = Gyx - YW*AWX
     +          +2.0D0*h02*AWXY
     +          +2.0D0*h0*(AWX-YW*AWXY)
      Gyy = Gyy + AW-YW*AWY
     +          -2.0D0*h02*AWYY
     +          +2.0D0*h0*YW*AWYY

c---
      if(Iopt.gt.1) then
c---

      cf = wn**3/(2.0D0*DW3)

      AWXYY= SX *(CHW*CHW+CHW*CX-2.0D0)*cf
      AWXXY=-SHW*( CX*CX +CHW*CX-2.0D0)*cf

      PWX = 2.0D0*(-AWX-2.0D0*h0*AWXY)
      PWY = 2.0D0*(-AWY+2.0D0*h0*AWYY)

      SXXX = AWX+YW*AWXY
     +              -2.0D0*h02*AWXYY
     +              +2.0D0*h0*YW*AWXYY
      SXYX = YW*AWYY
     +              -2.0D0*h02*AWXXY
     +              -2.0D0*h0*(AWYY-YW*AWXXY)
      SYXX = YW*AWYY
     +              +2.0D0*h02*AWXXY
     +              -2.0D0*h0*(AWYY+YW*AWXXY)
      SYYX = AWX-YW*AWXY
     +              -2.0D0*h02*AWXYY
     +              +2.0D0*h0*YW*AWXYY

      SXXY = 2.0D0*AWY+YW*AWYY
     +              +2.0D0*h0*(AWYY-YW*AWXXY)
     +              +2.0D0*h02*AWXXY
      SXYY = -AWX-YW*AWXY
     +              -2.0D0*h02*AWXYY
     +              +2.0D0*h0*(2.0D0*AWXY+YW*AWXYY)
      SYXY = -AWX-YW*AWXY
     +              +2.0D0*h02*AWXYY
     +              -2.0D0*h0*YW*AWXYY
      SYYY=       -YW*AWYY
     +              +2.0D0*h02*AWXXY
     +              +2.0D0*h0*(AWYY-YW*AWXXY)

      TXXX = - PWX + 2.0D0*SXXX
      TXXY =        SXXY+SYXX
      TYXX =   TXXY
      TYXY = - PWX + 2.0D0*SYXY

      TXYX = - PWY + 2.0D0*SXYX
      TXYY =        SXYY+SYYX
      TYYX =   TXYY
      TYYY = - PWY + 2.0D0*SYYY

      TXXX = TXXX + T1
      TXXY = TXXY + T2
      TYXX = TYXX + T2
      TYXY = TYXY + T3

      TXYX = TXYX + T2
      TXYY = TXYY + T3
      TYYX = TYYX + T3
      TYYY = TYYY + T4

      Px = Px+PWX
      Py = Py+PWY

c-----------
      end if
c-----------

c-----
c done
c-----

      return
      end
