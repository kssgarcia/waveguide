      subroutine slp_trgl3_sing
     +
     +   (x1,y1,z1
     +   ,x2,y2,z2
     +   ,x3,y3,z3
     +   ,GExx,GExy,GExz
     +   ,GEyx,GEyy,GEyz
     +   ,GEzx,GEzy,GEzz
     +   ,NGL
     +   ,Iflow
     +   )

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c---------------------------------------------
c Integrates the Green's function over a flat
c (three-node) triangle in local polar coordinates
c with origin at the singular point:
c
c (x1,y1,z1)
c
c SYMBOLS:
c -------
c
c asm: triangle area computed by numerical integration
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension zz(20),ww(20)

      common/var/wall

      common/zwl/zz,ww

      common/sgfr/a11,a12,a13,a21,a22,a23,a31,a32,a33
     +           ,b11,b12,b13,b21,b22,b23,b31,b32,b33
     +           ,ew,cell_vlm

      common/sgfi/Max1,Max2,Method
      common/NsNp/Ns,Np

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      piq = 0.25D0*pi
      pi4 = 4.0D0*pi

c------
c flags
c------

      Iopt = 1      ! for the Green's function call

c---
c compute surface metric: hs
c---

      dx = Dsqrt( (x2-x1)**2+(y2-y1)**2+(z2-z1)**2 )
      dy = Dsqrt( (x3-x1)**2+(y3-y1)**2+(z3-z1)**2 )

      vnx = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)
      vny = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)
      vnz = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)

      area = 0.5D0*sqrt(vnx*vnx + vny*vny + vnz*vnz)

      hs = 2.0D0*area   ! surface metric on a flat triangle

c---
c initialize
c---

      asm = 0.0D0    ! triangle area

      uxx = 0.0D0
      uxy = 0.0D0
      uxz = 0.0D0

      uyx = 0.0D0
      uyy = 0.0D0
      uyz = 0.0D0

      uzx = 0.0D0
      uzy = 0.0D0
      uzz = 0.0D0

c---
c apply the double quadrature
c---

      Do 1 i=1,NGL                  ! integration wrt phi

      ph    = piq*(1.0D0+zz(i))
      cph   = Dcos(ph)
      sph   = Dsin(ph)
      rmax  = 1.0D0/(cph+sph)
      rmaxh = 0.5D0*rmax

      sxx = 0.0D0
      sxy = 0.0D0
      sxz = 0.0D0

      syx = 0.0D0
      syy = 0.0D0
      syz = 0.0D0

      szx = 0.0D0
      szy = 0.0D0
      szz = 0.0D0

      bsm = 0.0D0     ! derivative of asm

          Do j=1,NGL               ! integration wrt r

          r  = rmaxh*(1.0D0+zz(j))
          xi = r*cph
          et = r*sph
          zt = 1.0D0-xi-et

          x = x1*zt + x2*xi + x3*et
          y = y1*zt + y2*xi + y3*et
          z = z1*zt + z2*xi + z3*et

c---
          if(Iflow.eq.1) then
c---

          call sgf_3d_fs
     +
     +       (Iopt
     +       ,x,y,z
     +       ,x1,y1,z1
     +       ,Gxx,Gxy,Gxz
     +       ,Gyx,Gyy,Gyz
     +       ,Gzx,Gzy,Gzz
     +       ,px,py,pz
     +       ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +       ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +       ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +       )

c---
          else if(Iflow.eq.2.or.Iflow.eq.6) then
c---

c         call sgf_3d_w    ! wall at x=wall
c    +
c    +       (Iopt
c    +       ,x,y,z
c    +       ,x1,y1,z1
c    +       ,wall
c    +       ,Gxx,Gxy,Gxz
c    +       ,Gyx,Gyy,Gyz
c    +       ,Gzx,Gzy,Gzz
c    +       ,px,py,pz
c    +       ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
c    +       ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
c    +       ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
c    +       )

          call sgf_3d_w    ! wall at y=wall
     +
     +       (Iopt
     +       ,y,z,x
     +       ,y1,z1,x1
     +       ,wall
     +       ,Gyy,Gyz,Gyx
     +       ,Gzy,Gzz,Gzx
     +       ,Gxy,Gxz,Gxx
     +       ,py,pz,px
     +       ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +       ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +       ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +       )

c---
          else if(Iflow.eq.3) then
c---

          call sgf_3d_3p
     +
     +        (Iopt
     +        ,x,y,z
     +        ,x1,y1,z1
     +        ,a11,a12,a13
     +        ,a21,a22,a23
     +        ,a31,a32,a33
     +        ,b11,b12,b13
     +        ,b21,b22,b23
     +        ,b31,b32,b33
     +        ,ew,cell_vlm
     +        ,Max1,Max2
     +        ,Gxx,Gxy,Gxz
     +        ,Gyx,Gyy,Gyz
     +        ,Gzx,Gzy,Gzz
     +        ,px,py,pz
     +        ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +        ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +        ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +        )

c---
          else if(Iflow.eq.7) then
c---

          call sgf_3d_2p_w    ! wall at y=wall
     +
     +      (Iopt
     +      ,y,z,x
     +      ,y1,z1,x1
     +      ,wall
     +      ,a11,a12,a21,a22
     +      ,Ns,Np
     +      ,Gyy,Gyz,Gyx
     +      ,Gzy,Gzz,Gzx
     +      ,Gxy,Gxz,Gxx
     +      ,py,pz,px
     +      ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +      ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +      ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +      )

c---
          else
c---
           write (6,*) " slp_trgl6_sing: invalid selection"
           stop
c---
          end if
c---

          cf = r*ww(j)

          bsm = bsm + cf

          sxx = sxx + Gxx*cf
          sxy = sxy + Gxy*cf
          sxz = sxz + Gxz*cf

          syx = syx + Gyx*cf
          syy = syy + Gyy*cf
          syz = syz + Gyz*cf

          szx = szx + Gzx*cf
          szy = szy + Gzy*cf
          szz = szz + Gzz*cf

      End Do

      cf = ww(i)*rmaxh

      asm = asm + bsm*cf

      uxx = uxx + sxx*cf
      uxy = uxy + sxy*cf
      uxz = uxz + sxz*cf
      uyx = uyx + syx*cf
      uyy = uyy + syy*cf
      uyz = uyz + syz*cf
      uzx = uzx + szx*cf
      uzy = uzy + szy*cf
      uzz = uzz + szz*cf

  1   Continue

c------------------------
c complete the quadrature
c------------------------

      cf = piq*hs

      asm  = asm*cf

      GExx = GExx + cf*uxx
      GExy = GExy + cf*uxy
      GExz = GExz + cf*uxz
      GEyx = GEyx + cf*uyx
      GEyy = GEyy + cf*uyy
      GEyz = GEyz + cf*uyz
      GEzx = GEzx + cf*uzx
      GEzy = GEzy + cf*uzy
      GEzz = GEzz + cf*uzz

c--------------------------
c  If all went well,
c  asm should be equal to "area"
c
c  write (6,100) i,area,asm
c--------------------------

c-----
c Done
c-----

 100  Format (1x,i3,2(f10.5))

      return
      end
