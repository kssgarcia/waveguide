      subroutine lgf_2d_crc
     + 
     +  (Iopt
     +  ,x,y
     +  ,x0,y0
     +  ,xc,yc
     +  ,a
     +  ,G
     +  ,Gx,Gy
     +  )

c=========================================
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c=========================================

c--------------------------------------------
c Neumman function of Laplace's equation
c in the exterior of a circle of radius ``a'' centered 
c at the point (xc, yc)
c
c Iopt =  1: compute only the Green's function
c      ne 1: compute the Green's function
c            and gradient
c
c see Pozrikidis (1997, p. 363)
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi

c---
c inverse point of the singularity
c---

      as  = a*a
      Dst = (x0-xc)**2+(y0-yc)**2
      xi  = xc + as*(x0-xc)/Dst
      yi  = yc + as*(y0-yc)/Dst

c----------------- 
c Neumman function
c---------------- 

      Dist0s = (x-x0)**2+(y-y0)**2     ! point sink
      Distis = (x-xi)**2+(y-yi)**2     ! image point source
      Distcs = (x-xc)**2+(y-yc)**2     ! center point sink

      G = log(Dist0s)+log(Distis)-log(Distcs)

      G = -0.5D0*G/pi2

      If(Iopt.eq.1) Go to 99

c-----
c compute the gradient
c----

      Gx = (x-x0)/Dist0s + (x-xi)/Distis
     +                   - (x-xc)/Distcs 
      Gy = (y-y0)/Dist0s + (y-yi)/Distis
     +                   - (y-yc)/Distcs 

      Gx = - Gx/pi2
      Gy = - Gy/pi2

c-----
c Done
c-----

  99  Continue

      Return
      End
