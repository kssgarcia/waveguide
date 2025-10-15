      subroutine Ilink 
     +
     +   (Nprtcl
     +   ,eps
     +   ,link
     +   )

c===========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c------------------------------------
c Compute the associaton matrix: link
c
c link(i,j) = 1 if particles i and j
c               are associated
c
c The matrix link is symmetric
c------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(25),y(25)
      Dimension link(25,25)

      common/points/x,y

c---
c associate with respect to distance
c---

      Do i=1,Nprtcl
       Do j=i+1,Nprtcl

        link(i,j) = 0

        dist = Dsqrt((x(i)-x(j))**2+(y(i)-y(j))**2)

        if(dist.lt.eps) link(i,j) = 1
        link(j,i) = link(i,j)

       End Do
      End Do

c-----
c Done
c-----

      return
      end
