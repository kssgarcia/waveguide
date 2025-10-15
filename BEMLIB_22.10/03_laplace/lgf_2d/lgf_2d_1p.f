      subroutine lgf_2d_1p
     +
     +   (Iopt
     +   ,x,y
     +   ,x0,y0
     +   ,RL
     +   ,G
     +   ,Gx,Gy
     +   )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c=======================================
c Singly periodic Green's function
c of Laplace's equation
c in two dimensions
c
c SYMBOLS:
c -------
c
c x, y :	coordinates of the field point
c x0,y0:	coordinates of the singular point
c
c RL:		period
c
c Iopt  = 1 produces G
c Iopt ne 1 produces G and its gradient (Gx,Gy)
c
c Notes:
c ----
c
c As y-y0 tends to infinity,
c
c G --> -(y-y0)/(2*RL)
c=======================================

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

      wn = pi2/RL       ! wave number

      A = wn*(y-y0)
      B = wn*(x-x0)
      C = Dcosh(A)-Dcos(B)

c-----------------
c Green's function
c----------------- 

      G = -log(2.0D0*C)/pi4

      if(Iopt.eq.1) Go to 99

c--------------------------
c Green's function gradient
c--------------------------

      cf = -0.5D0/(RL*C)

      Gx = cf *  Dsin(B)
      Gy = cf * Dsinh(A)

c-----
c Done
c-----

  99  Continue

      Return
      End
