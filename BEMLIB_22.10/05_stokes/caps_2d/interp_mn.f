      subroutine interp_mn (NSG,c,s,cm)

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------------
c Lagrange interpolation of concentration
c at mid-nodes from end-nodes
c
c Interpolation is done with respect to
c arc length (s)
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension c(0:900),cm(0:900),s(900)
      Dimension x(0:900)

c--------
c prepare
c--------

      NSG1 = NSG+1
      NSG2 = NSG+2

      x(0) = s(NSG)-s(NSG1)

      Do i=1,NSG1
        x(i) = s(i)
      End Do

      x(NSG2) = x(NSG1)+x(2)

      c(0)    = c(NSG)
      c(NSG1) = c(1)
      c(NSG2) = c(2)

c---
c launch
c---

      Do i=1,NSG

       x1 = x(i-1)
       x2 = x(i)
       x3 = x(i+1)
       x4 = x(i+2)

       y1 = c(i-1)
       y2 = c(i)
       y3 = c(i+1)
       y4 = c(i+2)

       xi = 0.5D0*(s(i)+s(i+1))

       cm(i) = y1*(xi-x2)*(xi-x3)*(xi-x4)/((x1-x2)*(x1-x3)*(x1-x4))
     +       + y2*(xi-x1)*(xi-x3)*(xi-x4)/((x2-x1)*(x2-x3)*(x2-x4))
     +       + y3*(xi-x1)*(xi-x2)*(xi-x4)/((x3-x1)*(x3-x2)*(x3-x4))
     +       + y4*(xi-x1)*(xi-x2)*(xi-x3)/((x4-x1)*(x4-x2)*(x4-x3))

      End Do

      cm(0)    = cm(NSG)
      cm(NSG1) = cm(1)
      cm(NSG2) = cm(2)

c-----
c done
c-----

      return
      end
