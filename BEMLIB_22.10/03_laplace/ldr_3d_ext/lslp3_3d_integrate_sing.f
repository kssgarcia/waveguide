      subroutine lslp3_3d_integrate_sing
     +
     +    (NGL
     +    ,x1,y1,z1
     +    ,x2,y2,z2
     +    ,x3,y3,z3
     +
     +    ,ztx1,zty1,ztz1
     +    ,ztx2,zty2,ztz2
     +    ,ztx3,zty3,ztz3
     +
     +    ,uxel,uyel,uzel
     +    ,area
     +    )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------------------------
c Compute the laplace single-layer potential
c for a vector density function over
c a linear (flat) triangle defined by three points 1-2-3
c
c Use local polar coordinates with origin at point 1
c
c SYMBOLS:
c -------
c
c asm: triangle area computed by numerical integration
c-------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension zz(20),ww(20)

      common/var1/Iflow,Ign,wall

      common/zwl/zz,ww
      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c------
c flags
c------

      Iopt_lgf = 1    ! need G only

c---
c compute triangle area and surface metric
c---

      dx = Dsqrt( (x2-x1)**2+(y2-y1)**2+(z2-z1)**2 )
      dy = Dsqrt( (x3-x1)**2+(y3-y1)**2+(z3-z1)**2 )

      vnx = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1)
      vny = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1)
      vnz = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)

      area = 0.5D0*Dsqrt( vnx**2+vny**2+vnz**2 )
      hs   = 2.0D0*area

      vnx = vnx/hs
      vny = vny/hs
      vnz = vnz/hs

c-----------
c initialize
c-----------

      asm = 0.0D0

      uux = 0.0D0
      uuy = 0.0D0
      uuz = 0.0D0

c---
c Double Gaussian quadrature
c---

      Do 1 i=1,NGL

      ph    = piq*(1.0D0+zz(i))
      cph   = Dcos(ph)
      sph   = Dsin(ph)
      rmax  = 1.0D0/(cph+sph)
      rmaxh = 0.5D0*rmax

      bsm = 0.0D0
      bix = 0.0D0
      biy = 0.0D0
      biz = 0.0D0

      Do 2 j=1,NGL

       r = rmaxh*(1.0D0+zz(j))

       xi = r*cph
       et = r*sph
       zt = 1.0D0-xi-et

       x = x1*zt + x2*xi + x3*et
       y = y1*zt + y2*xi + y3*et
       z = z1*zt + z2*xi + z3*et

       zetx = ztx1*zt + ztx2*xi + ztx3*et
       zety = zty1*zt + zty2*xi + zty3*et
       zetz = ztz1*zt + ztz2*xi + ztz3*et

c--------
       if(Iflow.eq.1) then
c--------

       call lgf_3d_fs 
     +
     +   (Iopt_lgf
     +   ,x,y,z
     +   ,x1,y1,z1
     +   ,G
     +   ,Gx,Gy,Gz
     +   )

c--------
       else if(Iflow.eq.2) then
c--------

      call lgf_3d_w
     +
     +   (Iopt_lgf
     +   ,Ign
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,wall
     +   ,G
     +   ,Gx,Gy,Gz
     +   )

c--------
      end if
c--------

       cf = r*ww(j)

       bsm = bsm+cf

       cf = cf*G
       bix = bix + zetx*cf
       biy = biy + zety*cf
       biz = biz + zetz*cf

  2   Continue

      cf = ww(i)*rmaxh

      asm = asm + bsm*cf

      uux = uux + bix*cf
      uuy = uuy + biy*cf
      uuz = uuz + biz*cf

  1   Continue

c---
c finish up the quadrature
c---

      cf = piq*hs

      asm = asm*cf

      uxel = uux*cf
      uyel = uuy*cf
      uzel = uuz*cf

c-----------------------------
c  if all went well,
c  asm should be equal to area
c
c     write (6,100) i,area,asm
c-----------------------------

 100  Format (1x,i3,2(f10.5))

      return
      end
