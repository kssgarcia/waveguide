      subroutine sgf_2d_1p 
     +
     +   (Iopt
     +   ,X,Y
     +   ,X0,Y0
     +   ,period
     +   ,Ising
     +   ,Gxx,Gxy
     +   ,Gyx,Gyy
     +   ,px,py
     +   ,Txxx,Txxy,Tyxx,Tyxy
     +   ,Txyx,Txyy,Tyyx,Tyyy
     +   )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------- 
c  Green's function of two-dimensional Stokes flow
c  associated with an array
c  of point forces deployed along the x axis
c
c  One of the point forces is located at 
c  the point: (X0,Y0)
c
c  The field point is located at: (X,Y)
c
c  SYMBOLS:
c  -------
c
c  period: distance between two point forces
c
c  Iopt = 1  Computes only G
c       ne 1  Computes G, p, T
c
c  If Ising = 1 will compute G-S,
c               where S is the Stokeslet
c               similarly for the stress and pressure
c------------------------------------------------- 

      Implicit Double Precision (a-h,o-z)
 
c---
c constants
c---

      pi  = 3.14159 265358D0
      pi2 = 2.0D0*pi

c----------------------------
c scale by the wave number
c----------------------------

      wn = pi2/period

      X  = X *wn   ! field point
      Y  = Y *wn
      X0 = X0*wn   ! singular point
      Y0 = Y0*wn

c----------
c compute G
c----------

      DX = X-X0
      DY = Y-Y0

      CHY = Dcosh(DY)
      SHY = Dsinh(DY)
      CX  = Dcos(DX)
      SX  = Dsin(DX)
      D   = CHY-CX 

      A  = 0.5D0*Dlog(2.0D0*D)
      AX = 0.5D0*SX /D
      AY = 0.5D0*SHY/D

      Gxx = -A-DY*AY   ! +1.0D0  ! optional
      Gxy =  DY*AX
      Gyx =  Gxy
      Gyy = -A+DY*AY

c---
c subtract off the Stokeslet
c if desired
c---

      If(Ising.eq.1) then

        DX2   = DX*DX
        DY2   = DY*DY
        DXY   = DX*DY
        DIST2 = DX2 + DY2
        SING  = 0.5D0*Dlog(DIST2)

        Gxx = Gxx + SING - DX2/DIST2
        Gxy = Gxy        - DXY/DIST2
        Gyx = Gyx        - DXY/DIST2 
        Gyy = Gyy + SING - DY2/DIST2

      End If

c---------------------------------
c  compute the stress and pressure
c  tensors
c---------------------------------

      if(Iopt.gt.1) then

      D2  = D*D
      AYY = 0.5D0*(1.0D0-CX*CHY)/D2
      AXY =-0.5D0*       SX*SHY /D2

      T1 = -2.0D0*(2.0D0*AX+DY*AXY)
      T2 = -2.0D0*(AY+DY*AYY)
      T3 =  2.0D0*DY*AXY
      T4 = -2.0D0*(AY-DY*AYY)

      px = 2.0D0*AX
      py = 2.0D0*AY

c---
c  subtract off the Stokeslet
c---

      If(Ising.eq.1) then

          cf  = 4.0D0/DIST2**2
          T1 = T1 + cf*DX*DX*DX
          T2 = T2 + cf*DX*DX*DY
          T3 = T3 + cf*DX*DY*DY
          T4 = T4 + cf*DY*DY*DY

          px = px - 2.0D0*DX/DIST2
          py = py - 2.0D0*DY/DIST2

      End If

c-------
c finish
c-------

      Txxx = T1*wn
      Txxy = T2*wn
      Tyxx = T2*wn
      Tyxy = T3*wn

      Txyx = T2*wn
      Txyy = T3*wn
      Tyyx = T3*wn
      Tyyy = T4*wn

      px = px*wn
      py = py*wn

      End If

c---
c reset
c---

      X  = X /wn
      Y  = Y /wn
      X0 = X0/wn
      Y0 = Y0/wn

c-----
c Done
c-----

      Return
      End
