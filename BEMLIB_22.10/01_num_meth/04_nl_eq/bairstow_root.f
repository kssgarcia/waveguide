      subroutine bairstow_root (a,b,c,x,n)

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c-----------------------------
c find the roots of a binomial
c-----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(20,2)

c---
c discriminant
c---

      disc = b**2-4.0D0*a*c
      a2   = 2.0D0*a

c-----
c quadratic formula
c-----

      if(disc.lt.0) then
        tmp  = Dsqrt(abs(disc))
        x(n,1)  = -b/a2
        x(n,2)  =  tmp/a2
        x(n-1,1)= -b/a2
        x(n-1,2)= -tmp/a2
      else
        tmp  = Dsqrt(disc)
        x(n,1)  = (-b+tmp)/a2
        x(n,2)  = 0.0D0
        x(n-1,1)= (-b-tmp)/a2
        x(n-1,2)= 0.0D0
      end if

c-----
c done
c-----

      return
      end
