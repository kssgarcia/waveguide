      subroutine surfact_fvm_etn
     +
     +  (Npts
     +  ,cel
     +  ,c
     +  )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------
c Compute node concentration from element
c concentration by averaging over neighbors
c------------------------------------------

      Implicit double precision (a-h,o-z)

      Dimension  p(1026,3)
      Dimension ne(1026,7)
      Dimension  c(1026)

      Dimension  cel(512)

      Dimension arel(512),xmom(512),ymom(512),zmom(512)

c--------------
c common blocks
c--------------

      common/points/p,ne

      common/geo1/arel
      common/geo9/xmom,ymom,zmom

c----------------
c Run over nodes
c---------------

      Do i=1,Npts

       collect = 0.0D0
       reduce  = 0.0D0

       Do j=1,ne(i,1)    ! run over adjacent elements

        k = ne(i,j+1)    ! element label

        xcnt = xmom(k)/arel(k)   ! element center
        ycnt = ymom(k)/arel(k)
        zcnt = zmom(k)/arel(k)

c---
c distance of ith node from element centroid
c---
        Dist = Dsqrt( (p(i,1)-xcnt)**2
     +               +(p(i,2)-ycnt)**2
     +               +(p(i,3)-zcnt)**2 )

        collect = collect + cel(k)/Dist
        reduce  = reduce  +  1.0D0/Dist

       End Do

       c(i) = collect/reduce

      End Do

c-----
c done
c-----

      return
      end
