        subroutine ttf
     +
     +   (n1
     +   ,x,y
     +   ,ux,uy
     +   ,ttfto
     +   ,ttftn
     +   )

c=============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=============================================

      Implicit Double Precision (a-h,o-z)

      Dimension x(100),y(100),ux(100),uy(100),um(100)
      Dimension ut(100),s(100)
      Dimension tnx(100),tny(100)

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c--------
c prepare
c--------

      n  = n1-1     ! number of intervals
      n2 = n+2

c------------
c wrap around
c------------

      x(n1) = x(1)
      y(n1) = y(1)

      ux(n1) = ux(1)
      uy(n1) = uy(1)

      x(n2) = x(2)
      y(n2) = y(2)

      ux(n2) = ux(2)
      uy(n2) = uy(2)

c---
c compute:  arc length
c---

      s(1) = 0.0

      Do i=2,n2
       ia = i-1
       s(i) = s(ia)+sqrt((x(i)-x(ia))**2
     +                  +(y(i)-y(ia))**2)
      End Do

c---
c compute:  tangential vector
c           by quadratic interpolation
c---

      Do i=2,n1

       ia = i-1
       i1 = i+1

       x0   = s(ia)-s(i)
       x1   = s(i1)-s(i)
       y0   = x(ia)-x(i)
       y1   = x(i1)-x(i)
       DxDs = (x0*y1/x1 - x1*y0/x0)/(x0-x1)
       y0   = y(ia)-y(i)
       y1   = y(i1)-y(i)
       DyDs = (x0*y1/x1 - x1*y0/x0)/(x0-x1)

       DmDs = sqrt(DxDs**2+DyDs**2)
       tnx(i) = DxDs/DmDs
       tny(i) = DyDs/DmDs

      End Do

      tnx(1) = tnx(n1)
      tny(1) = tny(n1)

c---
c compute:  velocity magnitude
c           tangential velocity
c---

      Do i=1,n1
       um(i) = sqrt(ux(i)**2+uy(i)**2)
       ut(i) = ux(i)*tnx(i)+uy(i)*tny(i)
      End Do

c---
c compute period by the trapezoidal rule
c---

      ttpto = 0.0
      ttptn = 0.0

      Do i=1,n
       i1 = i+1
       Ds = s(i1)-s(i)
       ttpto = ttpto + Ds*0.5*(1.0/um(i)+1.0/um(i1))
       ttptn = ttptn + Ds*0.5*(1.0/ut(i)+1.0/ut(i1))
      End Do

c---
c frequency
c----

      ttfto = pi4/ttpto
      ttftn = pi4/ttptn

c-----
c done
c-----

      return
      end
