      subroutine sslp_3d_integral_sing
     +
     +   (NGL
     +
     +   ,x1,y1,z1
     +   ,x2,y2,z2
     +   ,x3,y3,z3
     +
     +   ,fx1,fy1,fz1
     +   ,fx2,fy2,fz2
     +   ,fx3,fy3,fz3
     +
     +   ,uxel,uyel,uzel
     +   )

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c--------------------------------------------------------
c computes the single-layer potential over a linear
c (flat) triangle defined by three nodes 1-2-3
c
c Integration is done in local polar coordinates
c with origin at point 1
c
c SYMBOLS:
c -------
c
c asm: triangle area computed by numerical integration
c-------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension zz(20),ww(20)

      common/var/shrt,wall

      common/zwl/zz,ww
      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c---
c for triply periodic flow
c---

      common/aaaa/a11,a12,a13,a21,a22,a23,a31,a32,a33
      common/bbbb/b11,b12,b13,b21,b22,b23,b31,b32,b33
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

c------
c flags
c------

      Iopt_sgf = 1    ! only G is needed

c----------------------
c compute triangle area 
c and the surface metric
c----------------------

      vnx = (y2-y1)*(z3-z1) - (z2-z1)*(y3-y1)
      vny = (z2-z1)*(x3-x1) - (x2-x1)*(z3-z1)
      vnz = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)

      area = 0.5D0*Dsqrt(vnx*vnx +vny*vny+ vnz*vnz)
      hs   = 2.0D0*area

c--------------  ! for testing
c     uxx = 0.0
c     uxy = 0.0
c     uxz = 0.0
c     uyx = 0.0
c     uyy = 0.0
c     uyz = 0.0
c     uzx = 0.0
c     uzy = 0.0
c     uzz = 0.0
c--------------

      asm = 0.0D0

      uux = 0.0D0
      uuy = 0.0D0
      uuz = 0.0D0

      Do i=1,NGL

      ph    = piq*(1.0D0+zz(i))
      cph   = Dcos(ph)
      sph   = Dsin(ph)
      rmax  = 1.0D0/(cph+sph)
      rmaxh = 0.5D0*rmax

c--------------
c     sxx = 0.0
c     sxy = 0.0
c     sxz = 0.0
c     syx = 0.0
c     syy = 0.0
c     syz = 0.0
c     szx = 0.0
c     szy = 0.0
c     szz = 0.0
c--------------

      bsm = 0.0D0

      ssx = 0.0D0
      ssy = 0.0D0
      ssz = 0.0D0

      Do j=1,NGL

       r  = rmaxh*(1.0D0+zz(j))
       xi = r*cph
       et = r*sph
       zt = 1.0D0-xi-et

       x = x1*zt + x2*xi + x3*et
       y = y1*zt + y2*xi + y3*et
       z = z1*zt + z2*xi + z3*et

       dfx = fx1*zt + fx2*xi + fx3*et
       dfy = fy1*zt + fy2*xi + fy3*et
       dfz = fz1*zt + fz2*xi + fz3*et

       call sgf_3d_w 
     +
     +  (Iopt_sgf
     +  ,y,z,x
     +  ,y1,z1,x1
     +  ,wall
     +  ,gyy,gyz,gyx
     +  ,gzy,gzz,gzx
     +  ,gxy,gxz,gxx
     +  ,px,py,pz
     +  ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +  ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +  ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +  )

       cf = r*ww(j)

c---------------------------
c         sxx = sxx + gxx*cf
c         sxy = sxy + gxy*cf
c         sxz = sxz + gxz*cf
c         syx = syx + gyx*cf
c         syy = syy + gyy*cf
c         syz = syz + gyz*cf
c         szx = szx + gzx*cf
c         szy = szy + gzy*cf
c         szz = szz + gzz*cf
c---------------------------

       bsm = bsm + cf
       ssx = ssx + (dfx*gxx + dfy*gyx + dfz*gzx)*cf
       ssy = ssy + (dfx*gxy + dfy*gyy + dfz*gzy)*cf
       ssz = ssz + (dfx*gxz + dfy*gyz + dfz*gzz)*cf

      End Do

      cf = ww(i)*rmaxh

c-----------------------
c     uxx = uxx + sxx*cf
c     uxy = uxy + sxy*cf
c     uxz = uxz + sxz*cf
c     uyx = uyx + syx*cf
c     uyy = uyy + syy*cf
c     uyz = uyz + syz*cf
c     uzx = uzx + szx*cf
c     uzy = uzy + szy*cf
c     uzz = uzz + szz*cf
c-----------------------

      asm = asm + bsm*cf

      uux = uux + ssx*cf
      uuy = uuy + ssy*cf
      uuz = uuz + ssz*cf

      End Do

      cf = piq*hs

      asm = asm*cf

      uxel = uxel + uux*cf
      uyel = uyel + uuy*cf
      uzel = uzel + uuz*cf

c-----------------------------
c  if all went well,
c  asm should be equal to area
c
c     write (6,100) i,area,asm
c-----------------------------

 100  Format (1x,i3,2(f10.5))

      return
      end
