      subroutine lgf_2d_ww
     + 
     +  (Iopt
     +  ,Ign
     +  ,x,y
     +  ,x0,y0
     +  ,wall1,wall2
     +  ,G
     +  ,Gx,Gy
     +  )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c---------------------------------------------------
c Green's function of Laplace's equation in a domain
c confined between two parallel walls
c
c The first  wall is located at:  y = wall1
c The second wall is located at:  y = wall2
c
c Iopt =  1: compute only the Green's function
c      ne 1: compute  the Green's function and gradient
c
c (see Pozrikidis, 1997, p. 364)
c
c Set Ign = 1 to compute the green function
c     Ign = 2 to compute the Neumann function 
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.141592 65358D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c--------
c prepare
c--------

      h  = wall2-wall1
      h2 = 2.0D0*h

      wn = pi2/h2          ! wave number

      sing = 1.0D0                      ! Green
      if(Ign.eq.2) sing = -1.0D0        ! Neumann

c------------- 
c primary array
c------------- 

      dx = x-x0
      dy = y-y0

      A = wn*dx
      B = wn*dy
      C = Dcosh(A)-Dcos(B)

c---
c image system
c---
      
      y0i = 2.0D0*wall1-y0

      dxi = dx
      dyi = y-y0i

      Ai = wn*dxi
      Bi = wn*dyi
      Ci = Dcosh(Ai)-Dcos(Bi)

      G = ( -log(C) + sing * log(Ci) )/pi4

      if(Iopt.eq.1) Go to 99

c---------------------
c compute the gradient
c---------------------

      cf = -1.0D0/(2.0D0 * h2)

      Gx = cf * (sinh(A)/C - sing * sinh(Ai)/Ci)
      Gy = cf * (sin (B)/C - sing *  sin(Bi)/Ci)

c-----
c done
c-----

  99  Continue

      Return
      End
