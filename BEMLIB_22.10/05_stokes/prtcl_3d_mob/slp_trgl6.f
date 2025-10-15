      subroutine slp_trgl6
     +
     +  (x0,y0,z0
     +  ,k
     +  ,GExx,GExy,GExz
     +  ,GEyx,GEyy,GEyz
     +  ,GEzx,GEzy,GEzz
     +  ,mint
     +  ,Iflow
     +  )

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c----------------------------------------
c Integrates the Green's function over a
c non-singular triangle numbered: k
c
c SYMBOLS:
c -------
c
c mint: order of triangle quadrature
c
c GE_ij: integrated ij component over the element
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension p (1026,3)
      Dimension ne(1026,7)

      Dimension n(512,6),nbe (512,3)

      Dimension alpha(512),beta(512),gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/var/wall

      common/trq/xiq,etq,wq

      common/sgfr/a11,a12,a13,a21,a22,a23,a31,a32,a33
     +           ,b11,b12,b13,b21,b22,b23,b31,b32,b33
     +           ,ew,cell_vlm

      common/sgfi/Max1,Max2,Method
      common/NsNp/Ns,Np

c----------
c constants
c----------

      pi  = 3.14159 265358D0
      pi4 = 4.0D0*pi

c------
c flags
c------

      Ichoose = 2   ! glag for interp_p: need the metric hs
      Iopt    = 1   ! flag for the Green functions

c--------
c prepare
c--------

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
c vertices of the kth triangle
c---

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        call interp_p 
     +
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
     +  ,alpha(k),beta(k),gamma(k)
     +  ,xi,eta
     +  ,x,y,z
     +  ,DxDxi,DyDxi,DzDxi
     +  ,DxDet,DyDet,DzDet
     +  ,vnx,vny,vnz
     +  ,hxi,het,hs
     +  ,Ichoose
     +  )

c---
       if(Iflow.eq.1) then
c---

       call sgf_3d_fs 
     +
     +  (Iopt
     +  ,x,y,z
     +  ,x0,y0,z0
     +  ,Gxx,Gxy,Gxz
     +  ,Gyx,Gyy,Gyz
     +  ,Gzx,Gzy,Gzz
     +  ,px,py,pz
     +  ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +  ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +  ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +  )

c---
       else If(Iflow.eq.2) then
c---
c semi-infinite flow bounded by a wall at y=wall

c      call sgf_3d_w   ! wall at x=wall
c    +
c    +   (Iopt
c    +   ,x,y,z
c    +   ,x0,y0,z0
c    +   ,wall
c    +   ,Gxx,Gxy,Gxz
c    +   ,Gyx,Gyy,Gyz
c    +   ,Gzx,Gzy,Gzz
c    +   ,px,py,pz
c    +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
c    +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
c    +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
c    +   )

       call sgf_3d_w    ! wall at y=wall: NOTE DIFFERENT ARGUMENTS
     +
     +     (Iopt
     +     ,y,z,x
     +     ,y0,z0,x0
     +     ,wall
     +     ,Gyy,Gyz,Gyx
     +     ,Gzy,Gzz,Gzx
     +     ,gxy,gxz,gxx
     +     ,py,pz,px
     +     ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +     ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +     ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +     )

c---
       else if(Iflow.eq.3) then
c---

        call sgf_3d_3p
     +
     +    (Iopt
     +    ,x,y,z
     +    ,x0,y0,z0
     +    ,a11,a12,a13
     +    ,a21,a22,a23
     +    ,a31,a32,a33
     +    ,b11,b12,b13
     +    ,b21,b22,b23
     +    ,b31,b32,b33
     +    ,ew,cell_vlm
     +    ,Max1,Max2
     +    ,Gxx,Gxy,Gxz
     +    ,Gyx,Gyy,Gyz
     +    ,Gzx,Gzy,Gzz
     +    ,px,py,pz
     +    ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +    ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +    ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +    )

c---
       else if(Iflow.eq.7) then
c---

       call sgf_3d_2p_w    ! wall at y=wall
     +
     +  (Iopt
     +  ,y,z,x
     +  ,y0,z0,x0
     +  ,wall
     +  ,a11,a12,a21,a22
     +  ,Ns,Np
     +  ,Gyy,Gyz,Gyx
     +  ,Gzy,Gzz,Gzx
     +  ,gxy,gxz,gxx
     +  ,py,pz,px
     +  ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +  ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +  ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +  )

c---
       else
c---
         write (6,*) " intrg_trgl: invalid selection"
         stop
c---
       end if
c------------

       fc = 0.5D0*hs*wq(i)

       GExx = GExx + gxx*fc
       GExy = GExy + gxy*fc
       GExz = GExz + gxz*fc
       GEyx = GEyx + gyx*fc
       GEyy = GEyy + gyy*fc
       GEyz = GEyz + gyz*fc
       GEzx = GEzx + gzx*fc
       GEzy = GEzy + gzy*fc
       GEzz = GEzz + gzz*fc

      End Do

c-----
c Done
c-----

      return
      End
