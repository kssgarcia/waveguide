      subroutine interp_en (NSG,cm,s,c)

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c----------------------------------------
c Lagrange interpolation of concentration 
c at end-nodes from mid-nodes
c
c interpolation is done
c with respect to arc length (s)
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension c(0:200),cm(0:200),s(200)
      Dimension x(0:200)

c---
c prepare
c---

      NSG1 = NSG+1
      NSG2 = NSG+2

      x(0) = - 0.5D0*(s(NSG1)-s(NSG))

      Do i=1,NSG
        x(i) = 0.5D0*(s(i)+s(i+1))
      End Do

      x(NSG1) = s(NSG1)+x(1)
      x(NSG2) = s(NSG1)+x(2)

      cm(0)    = cm(NSG)
      cm(NSG1) = cm(1)
      cm(NSG2) = cm(2)

c---
c launch
c---

      Do i=2,NSG1

       x1 = x(i-2)
       x2 = x(i-1)
       x3 = x(i)
       x4 = x(i+1)

       y1 = cm(i-2)
       y2 = cm(i-1)
       y3 = cm(i)
       y4 = cm(i+1)

       xi = s(i)

       c(i) = y1*(xi-x2)*(xi-x3)*(xi-x4)/((x1-x2)*(x1-x3)*(x1-x4))
     +      + y2*(xi-x1)*(xi-x3)*(xi-x4)/((x2-x1)*(x2-x3)*(x2-x4))
     +      + y3*(xi-x1)*(xi-x2)*(xi-x4)/((x3-x1)*(x3-x2)*(x3-x4))
     +      + y4*(xi-x1)*(xi-x2)*(xi-x3)/((x4-x1)*(x4-x2)*(x4-x3))

      End Do

      c(0)    = c(NSG)
      c(1)    = c(NSG1)
      c(NSG2) = c(2)

c-----
c done
c-----

      return
      end
