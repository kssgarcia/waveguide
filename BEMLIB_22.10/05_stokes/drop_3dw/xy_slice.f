      subroutine xy_slice 
     +
     +  (npts
     +  ,cx,cy,cz
     +  ,nxy,jxy
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

c---------------------------
c identify nodes with z = cz
c---------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  p(1026,3)
      Dimension ne(1026,7)

      Dimension n(512,6),nbe(512,3)
      Dimension jxy(100)

      Logical switch

      common/points/p,ne
      common/elmnts/n,nbe

c---
c introduce an integer array jxy that contains
c the labels of all points in the z = cz plane
c---

      nxy = 0   ! counter

      Do i=1,npts
        if(abs(p(i,3)-cz).lt.0.000001) then
          nxy = nxy+1
          jxy(nxy) = i
        end if
      End Do

c---
c sort marker points
c according to meridional angle
c using bubble-sort
c---

      switch = .true.

      Do while (switch)

        switch = .false.

        Do i=1,nxy-1
          i1 = i+1
          phii  = atan2(p(jxy(i ),1)-cx,
     +                  p(jxy(i ),2)-cy)
          phip1 = atan2(p(jxy(i1),1)-cx,
     +                  p(jxy(i1),2)-cy)

          if(phii.gt.phip1) then
             jhold = jxy(i)
             jxy(i)  = jxy(i1)
             jxy(i1) = jhold
             switch = .true.
            end if
         End Do

      End Do

c---
c close loop
c---

      nxy = nxy+1
      jxy(nxy) = jxy(1)

c-----
c done
c-----

      return
      end
