      subroutine lgf_3d_sph
     +
     +   (Iopt
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,xc,yc,zc
     +   ,a
     +   ,G
     +   ,Gx,Gy,Gz
     +   )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c--------------------------------------------
c Neumman function for Laplace's equation
c in a domain bounded internally by a sphere
c
c SYMBOLS:
c -------
c
c a: sphere radius
c
c xc,yc,zc: sphere center
c
c Iopt =  1 compute only the Green's function
c      ne 1 compute the Green's function and gradient
c
c NN: intervals for numerical integration
c     of the image sink distribution
c
c see Pozrikidis (1997, p. 335)
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (NN=50)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi4 = 4.0D0*pi

c-------------------- 
c primary singularity
c------------------- 

      dx = x-x0
      dy = y-y0
      dz = z-z0

      r = sqrt(dx**2+dy**2+dz**2)

      G = 1.0D0/(pi4*r)

c-------------------------------
c singularity at the image point
c-------------------------------

      as   = a**2
      Dists = (x0-xc)**2+(y0-yc)**2+(z0-zc)**2
      cf   = as/Dists
      x0i  = xc + cf * (xs-xc)
      y0i  = yc + cf * (ys-yc)
      z0i  = zc + cf * (zs-zc)

      dxi = x-x0i
      dyi = y-y0i
      dzi = z-z0i

      ri = Dsqrt(dxi**2+dyi**2+dzi**2)

      RR = Dsqrt(Dists)
      
      fc = a/(pi4*RR)

      G = G + fc * 1.0D0/(pi4*ri)

c---
c distribution of image point sinks
c---

      sum = 0.0D0       ! distribution

      Do i=1,NN     ! integration by the trapezoidal rule
        ff = (i-0.5D0)/NN
        xx = xc + (x0i-xc)*ff
        yy = yc + (y0i-yc)*ff
        zz = zc + (z0i-zc)*ff
        Dist = sqrt((x-xx)**2+(y-yy)**2+(z-zz)**2)
        sum = sum + 1.0D0/Dist
      End Do

      G = G - fc*sum/NN


c---------------------
c compute the gradient
c---------------------

      if(Iopt.gt.1) then

      den  = r**3  *pi4
      deni = ri**3 *pi4

      Gx = - Dx/den - fc* Dxi/deni
      Gy = - Dy/den - fc*Dyi/deni
      Gz = - Dz/den - fc*Dzi/deni

      sumx = 0.0D0      ! distribution
      sumy = 0.0D0
      sumz = 0.0D0

      Do i=1,NN     ! integration by the trapezoidal rule

        ff = (i-0.5D0)/NN
        xx = xc + (x0i-xc)*ff
        yy = yc + (y0i-yc)*ff
        zz = zc + (z0i-zc)*ff
        Distc = sqrt((x-xx)**2+(y-yy)**2+(z-zz)**2)**3
        sumx = sumx + (x-xx)/Distc
        sumy = sumy + (y-yy)/Distc
        sumz = sumz + (z-zz)/Distc

      End Do

      Gx = Gx + fc*sumx/NN
      Gy = Gy + fc*sumy/NN
      Gz = Gz + fc*sumz/NN

      End If

c-----
c done
c-----

      Return
      End
