      subroutine lgf_2d_fs 
     +
     +    (Iopt
     +    ,x,y
     +    ,x0,y0
     +    ,G
     +    ,Gx,Gy
     +    )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c----------------------------------------
c Free-space Green's function of Laplace's
c equation in two dimensions:
c
c  G = -1/(2*pi) * lnr
c
c  Iopt =  1 computes only the Green's function
c       ne 1 computes the Green's function
c                     and its gradient
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c constants
c---

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi

c--- 
c launch
c--- 

      dx = x-x0
      dy = y-y0

      rs = dx*dx+dy*dy

      G  = - 0.5D0*log(rs)/pi2

      if(Iopt.gt.1) then

       den = rs*pi2
       Gx = - dx/den
       Gy = - dy/den

      end if

c-----
c Done
c-----

      return
      end
