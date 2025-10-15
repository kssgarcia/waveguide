      subroutine lgf_2d_fs 
     +
     +  (Iopt
     +  ,x,y
     +  ,x0,y0
     +  ,G
     +  ,Gx,Gy
     +  )

c=========================================
c FDLIB, CFDLAB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=========================================

c-------------------------------------------
c Free-space Green's function:
c
c  G = -(1/2*pi) * lnr
c
c Iopt = 1: compute only the Green's function
c     ne 1: compute the Green's function
c           and the gradient
c
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

c-----------
c  constants
c-----------

      pi  = 3.14159 265358D0
      pi2 = 2.0D0*pi

c-----------------
c Green's function
c-----------------

      dx = x-x0
      dy = y-y0

      rs = dx**2+dy**2

      G  = - 0.5D0*log(rs)/pi2

      If(Iopt.eq.1) Go to 99

c---------
c Gradient
c---------

      den = rs*pi2

      Gx = - dx/den
      Gy = - dy/den

c-----
c Done
c-----

  99  Continue

      Return
      End
