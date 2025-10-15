      subroutine sgf_2d_1p_w 
     +
     +              (Iselect
     +              ,Xfl,Yfl
     +              ,Xpf,Ypf
     +              ,wall
     +              ,period
     +              ,Ising
     +              ,Gxx,Gxy
     +              ,Gyx,Gyy
     +              ,px,py
     +              ,Txxx,Txxy,Tyxx,Tyxy
     +              ,Txyx,Txyy,Tyyx,Tyyy
     +              )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------- 
c  Green's function of two-dimensional Stokes flow
c  associated with an array array of point forces 
c  deployed along the x axis above a wall
c
c  The wall is located at y = wall
c
c  One of the point forces is located at (Xpf,Ypf)
c
c  The field point is located at (Xfl,Yfl)
c
c  SYMBOLS:
c  -------
c
c  period:  Distance between the point forces
c
c  Iselect = 1 computes only G
c           ne computes G, p, T
c
c  Ising = 1 computes G-S where S is the Stokeslet
c            similarly for the stress and pressure 
c
c------------------------------------------------- 

      Implicit Double Precision (a-h,o-z)
 
c---
c constants
c---

      pi  = 3.14159 265358
      pi2 = 2.0*pi

c-----------------------------------------------------
c Change length scale to make the period equal to 2*pi
c-----------------------------------------------------

      scale = pi2/period

      X  = Xfl * scale   ! field point
      Y  = Yfl * scale

      X0 = Xpf * scale   ! singular point
      Y0 = Ypf * scale

      WP = wall * scale

c-------------------------------
c periodic array of point forces
c-------------------------------

      DX = X-X0
      DY = Y-Y0

      CHY = COSH(DY)
      SHY = SINH(DY)
      CX  = COS(DX)
      SX  = SIN(DX)
      D   = CHY-CX 

      A   = 0.5*LOG(2.0*D)
      AX  = 0.5*SX /D
      AY  = 0.5*SHY/D

      Gxx = -A-DY*AY+1.0
      Gxy =  DY*AX
      Gyx =  Gxy
      Gyy = -A+DY*AY

c---
c subtract off the Stokeslet
c---

      If(Ising.eq.1) then

        DX2   = DX**2
        DY2   = DY**2
        DXY   = DX*DY
        DIST2 = DX2+DY2
        SING  = 0.5*LOG(DIST2)

        Gxx = Gxx + SING - DX2/DIST2
        Gxy = Gxy        - DXY/DIST2
        Gyx = Gyx        - DXY/DIST2 
        Gyy = Gyy + SING - DY2/DIST2

      End If

      If(Iselect.eq.1) Go to 98

c----------------------------
c compute stress and pressure
c----------------------------

      D2  = D**2
      AYY = 0.5*(1.0-CX*CHY)/D2
      AXY =-0.5*     SX*SHY /D2

      T1 = -2.0*(2.0*AX+DY*AXY)
      T2 = -2.0*(AY+DY*AYY)
      T3 =  2.0*DY*AXY
      T4 = -2.0*(AY-DY*AYY)

      PX = 2.0*AX
      PY = 2.0*AY

c----
c subtract off the Stokeslet
c----

      If(Ising.eq.1) then

          DIST4 = DIST2**2
          cf    = 4.0/DIST4 

          T1 = T1 + cf*DX**3  
          T2 = T2 + cf*DX**2*DY
          T3 = T3 + cf*DX   *DY**2
          T4 = T4 + cf      *DY**3 

          PX = PX - 2.0*DX/DIST2
          PY = PY - 2.0*DY/DIST2

        End If

  98    Continue

c--------------------
c Image singularities
c--------------------

      YW = Y+Y0-2.0*WP

      CHW = cosh(YW)
      SHW = sinh(YW)
      DW  = CHW-CX
      DW2 = DW**2
      DW3 = DW2*DW

      AW   =  0.5*LOG(2.0*DW)
      AWX  =  0.5*SX /DW
      AWY  =  0.5*SHW/DW
      AWYY =  0.5*(1.0-CX*CHW)/DW2
      AWXY = -0.5*     SX*SHW /DW2

      RL  = Y0-WP
      RL2 = RL**2

      Gxx = Gxx + AW-1.0+YW*AWY
     +                        +2.0*RL*YW*AWYY
     +                                            -2.0*RL2*AWYY
      Gxy = Gxy - YW*AWX
     +                        +2.0*RL*(AWX+YW*AWXY)
     +                                            -2.0*RL2*AWXY
      Gyx = Gyx - YW*AWX
     +                        +2.0*RL*(AWX-YW*AWXY)
     +                                            +2.0*RL2*AWXY
      Gyy = Gyy + AW-YW*AWY
     +                        +2.0*RL*YW*AWYY
     +                                            -2.0*RL2*AWYY

c---
      If(Iselect.eq.1) Go to 99
c---

      cf = 1/(2.0*DW3)

      AWXYY= SX *(CHW**2+CHW*CX-2.0)*cf
      AWXXY=-SHW*( CX**2+CHW*CX-2.0)*cf

      PWX = 2.0*(-AWX-2.0*RL*AWXY)
      PWY = 2.0*(-AWY+2.0*RL*AWYY)

      SXXX = AWX+YW*AWXY
     +              +2.0*RL*YW*AWXYY
     +                                    -2.0*RL2*AWXYY
      SXYX = YW*AWYY
     +              -2.0*RL*(AWYY-YW*AWXXY)
     +                                    -2.0*RL2*AWXXY
      SYXX = YW*AWYY
     +              -2.0*RL*(AWYY+YW*AWXXY)
     +                                    +2.0*RL2*AWXXY
      SYYX = AWX-YW*AWXY
     +              +2.0*RL*YW*AWXYY
     +                                    -2.0*RL2*AWXYY

      SXXY = 2.0*AWY+YW*AWYY
     +              +2.0*RL*(AWYY-YW*AWXXY)
     +                                    +2.0*RL2*AWXXY
      SXYY = -AWX-YW*AWXY
     +              +2.0*RL*(2.0*AWXY+YW*AWXYY)
     +                                    -2.0*RL2*AWXYY
      SYXY = -AWX-YW*AWXY
     +              -2.0*RL*YW*AWXYY
     +                                    +2.0*RL2*AWXYY
      SYYY=       -YW*AWYY
     +                 +2.0*RL*(AWYY-YW*AWXXY)
     +                                    +2.0*RL2*AWXXY

      TXXX = - PWX + 2.0*SXXX
      TXXY =        SXXY+SYXX
      TYXX =   TXXY
      TYXY = - PWX + 2.0*SYXY

      TXYX = - PWY + 2.0*SXYX
      TXYY =        SXYY+SYYX
      TYYX =   TXYY
      TYYY = - PWY + 2.0*SYYY

      TXXX = TXXX + T1
      TXXY = TXXY + T2
      TYXX = TYXX + T2
      TYXY = TYXY + T3

      TXYX = TXYX + T2
      TXYY = TXYY + T3
      TYYX = TYYX + T3
      TYYY = TYYY + T4

      PX = PX+PWX
      PY = PY+PWY

c-----------
c scale back
c-----------

      TXXX = TXXX*scale
      TXXY = TXXY*scale
      TYXX = TYXX*scale
      TYXY = TYXY*scale

      TXYX = TXYX*scale
      TXYY = TXYY*scale
      TYYX = TYYX*scale
      TYYY = TYYY*scale

      PX = PX*scale
      PY = PY*scale

c---
c Done
c---

  99  Continue

      Return
      End
