      subroutine lgf_2d_www
     + 
     +  (Iopt
     +  ,x,y
     +  ,x0,y0
     +  ,wall1,wall2,wall3
     +  ,G
     +  ,Gx,Gy
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c Neumman function of Laplace's equation in a 
c semi-infinite strip confined by three walls:
c
c a plane wall located at x = wall1
c a plane wall located at x = wall2
c a plane wall located at y = wall3
c
c Iopt =  1: compute only the Green's function 
c      ne 1: compute  the Green's function and gradient
c
c For formulae see Pozrikidis (1997, p. 364)
c--------------------------------------------

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

      h = wall2-wall1
      h2 = 2.0D0*h
      wn = pi2/h2          ! wave number

c--------------------------- 
c singularity and image system
c with respect to wall1
c--------------------------- 

      dx  = x-x0
      dxi = x-2.0D0*wall1+x0
      A   = wn*dx
      Ai  = wn*dxi

      dy = y-y0
      B  = wn*dy
      C  = cosh(B)
      D  = C-cos(A)
      Di = C-cos(Ai)

      G = -log(D*Di)/pi4

c-----------------------------
c Images with respect to wall3
c-----------------------------

      dy  = y-2.0*wall3+y0
      BB  = wn*dy
      CC  = cosh(BB)
      DD  = CC-cos(A)
      DDi = CC-cos(Ai)

      G = G-log(DD*DDi)/pi4

      If(Iopt.eq.1) Go to 99

c---------------------
c compute the gradient
c---------------------

      cf = -1.0/(2.0D0*h2)

      tmp  = 1.0D0/D +1.0D0/DD
      tmpi = 1.0D0/Di+1.0D0/DDi

      tlp  = 1.0D0/D +1.0D0/Di
      tlpi = 1.0D0/DD+1.0D0/DDi

      Gx = cf * (sin (A)*tmp + sin (Ai)*tmpi)
      Gy = cf * (sinh(B)*tlp + sinh(BB)*tlpi)

c-----
c Done
c-----

  99  Continue

      Return
      End
