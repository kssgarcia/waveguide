      subroutine lgf_2d_1p_w
     +
     +  (Iopt
     +  ,Ign
     +  ,x,y
     +  ,x0,y0
     +  ,RL
     +  ,wall
     +  ,G
     +  ,Gx,Gy
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------
c Periodic Green's function of Laplace's equation
c in a semi-infinite domain bounded by a plane wall.
c located at y = wall
c
c SYMBOLS:
c -------
c
c x, y :  coordinates of the field point
c x0,y0:  coordinates of the singular point
c RL:     period
c
c Iopt  = 1 produces G
c Iopt ne 1 produces G and its gradient (Gx,Gy)
c
c Ign = 1 produces the Green's function
c Ign = 2 produces the Neumann function
c
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c--------
c prepare
c--------

      sign = 1.0D0
      If(Ign.eq.2) sign = -1.0D0

      wn = pi2/RL       ! wave number

      xh = wn*(x-x0)
      yh = wn*(y-y0)

      C = Dcosh(yh)-Dcos(xh)

c-----------------
c Green's function
c---------------- 

      G = -log(2.0D0*C)/pi4

c--- image singularity:
c    -----------------

      x0i = x0
      y0i = 2.0D0*wall - y0

      xhi = xh
      yhi = wn*(y-y0i)

      Ci = cosh(yhi)-cos(xhi)

      G = G + sign*Dlog(2.0D0*Ci)/pi4

c--------------------------
c Green's function gradient
c--------------------------

      If(Iopt.eq.1) Go to 99

      cf  = -0.5D0/(RL*C)

      cfi = -0.5D0/(RL*Ci)

      Gx = cf* sin(xh) - sign*cfi* Dsin(xhi)

      Gy = cf*sinh(yh) - sign*cfi*Dsinh(yhi)

c-----
c Done
c-----

  99  Continue

      Return
      End
