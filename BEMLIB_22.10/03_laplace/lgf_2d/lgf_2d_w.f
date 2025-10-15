      subroutine lgf_2d_w
     +
     +  (Iopt
     +  ,Ign
     +  ,x,y
     +  ,x0,y0
     +  ,wall
     +  ,G
     +  ,Gx,Gy
     +  )

c==========================================
c FDLIB, CFDLAB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c==========================================

c---------------------------------------
c Green and Neumann functions of Laplace's equation
c in a semi-infinite domain bounded by a plane wall
c located at y = wall
c
c  The Green function is given by:
c
c  G = -1/(2*pi)*ln(r) + 1/(2*pi)*ln(r_image)
c
c  The Neumann function is given by:
c
c  G = -1/(2*pi)*ln(r) - 1/(2*pi)*ln(r_image)
c
c  LEGEND:
c  ------
c
c  Iopt =  1: compute only the Green's function
c       ne 1: compute the Green's function
c             and its gradient
c
c  Ign =  1: for the Green function
c         2: for the Neumann function
c
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

c--------
c prepare
c--------

      sign = 1.0D0
      If(Ign.eq.2) sign = -1.0D0

c-------------------- 
c primary singularity
c------------------- 

      dx = x-x0
      dy = y-y0

      rs = dx**2+dy**2

      G = - 0.5D0*Dlog(rs)/pi2

c--- image singularity:
c    -----------------

      x0i = x0
      y0i = 2.0D0*wall- y0

      dxi = dx
      dyi = y-y0i
      rsi = dxi**2+dyi**2

      G = G + 0.5D0 * sign*Dlog(rsi)/pi2

      If(Iopt.eq.1) Go to 99

c---------------------
c compute the gradient
c---------------------

      den  = rs *pi2
      deni = rsi*pi2

c--- image singularity:
c    -----------------

      Gx = - Dx/den + sign * Dxi/deni
      Gy = - Dy/den + sign * Dyi/deni

c-----
c Done
c-----

  99  Continue

      Return
      End
