      subroutine slp_trgl6_sing
     +
     +   (x0,y0,z0
     +   ,k
     +   ,GExx,GExy,GExz
     +   ,GEyx,GEyy,GEyz
     +   ,GEzx,GEzy,GEzz
     +   ,NGL
     +   ,Iflow
     +   )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c----------------------------------------
c Integrate the Green's function over the kth
c singular quadratic triangle
c
c This is done by breaking up the singular
c triangle into six flat triangles, and
c then integrating individually over the
c flat triangles in local poral coordinates
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension p (1026,3)
      Dimension ne(1026,7)
      Dimension n (512,6),nbe (512,3)

      common/points/p,ne
      common/elmnts/n,nbe

      common/var/wall

c---
c launching
c---

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      GExx = 0.0D0
      GExy = 0.0D0
      GExz = 0.0D0
      GEyx = 0.0D0
      GEyy = 0.0D0
      GEyz = 0.0D0
      GEzx = 0.0D0
      GEzy = 0.0D0
      GEzz = 0.0D0

c---
c Integrate over six flat triangles
c---

      call slp_trgl3_sing
     +
     +  (x0,y0,z0
     +  ,p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,GExx,GExy,GExz
     +  ,GEyx,GEyy,GEyz
     +  ,GEzx,GEzy,GEzz
     +  ,NGL,Iflow
     +  )

      call slp_trgl3_sing
     +
     +  (x0,y0,z0
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,GExx,GExy,GExz
     +  ,GEyx,GEyy,GEyz
     +  ,GEzx,GEzy,GEzz
     +  ,NGL,Iflow
     +  )

      call slp_trgl3_sing
     +
     +  (x0,y0,z0
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,GExx,GExy,GExz
     +  ,GEyx,GEyy,GEyz
     +  ,GEzx,GEzy,GEzz
     +  ,NGL,Iflow
     +  )

      call slp_trgl3_sing
     +
     +   (x0,y0,z0
     +   ,p(i5,1),p(i5,2),p(i5,3)
     +   ,p(i3,1),p(i3,2),p(i3,3)
     +   ,GExx,GExy,GExz
     +   ,GEyx,GEyy,GEyz
     +   ,GEzx,GEzy,GEzz
     +   ,NGL,Iflow
     +   )

      call slp_trgl3_sing
     +
     +   (x0,y0,z0
     +   ,p(i3,1),p(i3,2),p(i3,3)
     +   ,p(i6,1),p(i6,2),p(i6,3)
     +   ,GExx,GExy,GExz
     +   ,GEyx,GEyy,GEyz
     +   ,GEzx,GEzy,GEzz
     +   ,NGL,Iflow
     +   )

      call slp_trgl3_sing
     +
     +   (x0,y0,z0
     +   ,p(i6,1),p(i6,2),p(i6,3)
     +   ,p(i1,1),p(i1,2),p(i1,3)
     +   ,GExx,GExy,GExz
     +   ,GEyx,GEyy,GEyz
     +   ,GEzx,GEzy,GEzz
     +   ,NGL,Iflow
     +   )

c-----
c done
c-----

       return
       end
