      subroutine lgf_3d_fs 
     +
     +   (Iopt
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,G
     +   ,Gx,Gy,Gz
     +   )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c----------------------------------------
c Free-space Green's function of Laplace's equation:
c
c  G(x, x0) = 1/(4*pi*r)
c 
c  where: r = |x-x0|
c
c SYMBOLS:
c --------
c
c  Iopt =  1 compute the Green's function G
c       ne 1 compute the Green's function G
c            and its gradient (Gx, Gy, Gz)
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi

c------ 
c launch
c------ 

      Dx = x-x0
      Dy = y-y0
      Dz = z-z0

      r = Dsqrt(Dx*Dx + Dy*Dy + Dz*Dz)

      G = 1.0D0/(pi4*r)


c---------
c gradient
c---------

      if(iopt.gt.1) then

      den = pi4*r*r*r

      Gx = - Dx/den
      Gy = - Dy/den
      Gz = - Dz/den

      end if

c-----
c done
c-----

      return
      end
