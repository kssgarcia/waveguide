      subroutine cubic
     +
     +  (a,b,c
     +  ,D
     +  ,x1,x2,x3
     +  ,x23r,x23i
     +  )

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c================================
c Roots of the cubic equation:
c
c  x**3 + a*x**2 + b*x + c = 0
c
c where a,b,c are real coefficients.
c
c The roots are computed analytically
c using Cardano's analytical formulas
c
c CASE A
c
c three real roots: x1, x2, x3
c
c CASE B
c
c one real root: x1
c two complex conjugate roots: x23r +- x23i
c================================

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      oot   = 1.0D0/3.0D0
      opf   = 1.5D0
      three = 3.0D0
      srth  = sqrt(three)

c--------
c prepare
c--------

      p = (3.0D0*b-a*a)/3.0D0
      q = c+2.0D0*a*a*a/27.0D0-a*b/3.0D0

      D = (p/3.0D0)**3+(q/2.0D0)**2        ! discriminant

c--------------------------------------
c one real root: x1
c two complex conjugate roots: x23r +- x23i
c--------------------------------------

      if(D.ge.0) then 

       srd  = sqrt(D)
       tmp  =  -0.5D0*q+srd
       u    =  abs(tmp)**oot

       if(tmp.lt.0) u = -u

       tmp  =  -0.5D0*q-srd
       v    =  abs(tmp)**oot

       if(tmp.lt.0) v = -v

       x1   = -a/3.0D0+u+v

       x23r = -a/3.0D0-0.5D0*(u+v)
       x23i =   srth*0.5D0*(u-v)

c-----------------
c three real roots: x1, x2, x3
c-----------------

      else       ! three real roots

       cosphi = -0.5D0*q/(abs(p)/3.0D0)**opf
       phi    = acos(cosphi)

       cf = 2.0*sqrt(abs(p)/3.0D0)

       x1 = -a/3.0D0+cf*cos(phi/3.0D0)
       x2 = -a/3.0D0-cf*cos((phi-pi)/3.0D0)
       x3 = -a/3.0D0-cf*cos((phi+pi)/3.0D0)

c-----------
      end if
c-----------

c---
c Done
c---

      return
      end
