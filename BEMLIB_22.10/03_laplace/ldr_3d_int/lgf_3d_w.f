      subroutine lgf_3d_w
     +
     +   (Iopt
     +   ,Ign
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,wall
     +   ,G
     +   ,Gx,Gy,Gz
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

c---------------------------------------
c Green's function of Laplace's equation
c in a semi-infinite domain
c bounded by a plane wall located at x = wall
c
c  G = 1/(4*pi*r) (+-)  1/(4*pi*r_im) 
c
c Iopt =  1 computes only the Green's function
c      ne 1 computes the Green's function
c           and the gradient
c
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi

c--------
c prepare
c--------

      sign = 1.0D0
      If(Ign.eq.2) sign = -1.0D0

c--------------------
c primary singularity
c-------------------- 

      dx = x-x0
      dy = y-y0
      dz = z-z0

      r = Dsqrt(dx**2+dy**2+dz**2)

      G = 1.0D0/(pi4*r)

c------------------
c image singularity
c------------------

      x0i = 2.0D0*wall - x0
      y0i = y0
      z0i = z0

      dxi = x-x0i
      dyi = dy
      dzi = dz

      ri = Dsqrt(dxi**2+dyi**2+dzi**2)

      G = G - sign/(pi4*ri)

      If(Iopt.eq.1) Go to 99

c---------------------
c compute the gradient
c---------------------

      den  = r**3  *pi4
      deni = ri**3 *pi4

      Gx = - Dx/den + sign * Dxi/deni
      Gy = - Dy/den + sign * Dyi/deni
      Gz = - Dz/den + sign * Dzi/deni

c-----
c Done
c-----

  99  Continue

      Return
      End
