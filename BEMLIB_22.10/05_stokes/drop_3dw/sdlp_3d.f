      subroutine sdlp_3d
     +
     +  (npts
     +  ,nelm
     +  ,mint
     +  ,u
     +  ,dlp
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c-------------------------------------------------
c Compute the double-layer potential at the nodes
c
c Desingularization is done by use of an 
c integral identity
c
c rbase : radius of the circular base
c xid :   quadrature for a disk with mdsk points
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)
      Dimension    u(1026,3)

      Dimension nvel(1026)

      Dimension  dlp(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xid(30),etd(30),wd(30)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/var/shrt,wall
      common/veloc1/nvelt,nvel

      common/base/rbase
      common/dski/mdsk
      common/dskr/xid,etd,wd

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c--------
c prepare
c--------

      rbases = rbase**2  ! circular drop base radius squared

c--------------
c for debugging
c--------------

c     Do i=1,npts
c       u(i,1) = 0.0D0
c       u(i,2) = 0.0D0
c       u(i,3) = 0.0D0
c     End do
c
c---------------------

c----------------
c loop over nodes
c----------------

      Do 1 node=1,nvelt

      i = nvel(node)

      x0 = p(i,1)
      y0 = p(i,2)
      z0 = p(i,3)

      u0 = u(i,1)
      v0 = u(i,2)
      w0 = u(i,3)

c----------------------------------------
c     write (6,*) " Enter x0,y0,z0"      ! to debug
c     read  (5,*) x0,y0,z0               ! to debug
c     write (6,*) " Enter u0,v0,w0"      ! to debug
c     read  (5,*) u0,v0,w0               ! to debug
c----------------------------------------

c---
c compute the dlp over the interface
c---

      us = 0.0D0
      vs = 0.0D0
      ws = 0.0D0

      Do k=1,nelm

        call sdlp_3d_integral
     +
     +     (x0,y0,z0
     +     ,u0,v0,w0
     +     ,k
     +     ,mint
     +     ,u
     +     ,uxel,uyel,uzel
     +     )

        us = us + uxel
        vs = vs + uyel
        ws = ws + uzel

c       write (6,100) k,uxel,uyel,uzel

      End Do

c---
c Compute the dlp over the base
c---

      Iopt_sgf = 2

      sumxx = 0.0D0
      sumxy = 0.0D0
      sumxz = 0.0D0

      sumyx = 0.0D0
      sumyy = 0.0D0
      sumyz = 0.0D0

      sumzx = 0.0D0
      sumzy = 0.0D0
      sumzz = 0.0D0

      Do j=1,mdsk

      x = rbase * xid(j)
      y = 0.0D0
      z = rbase * etd(j)

      call sgf_3d_w
     +
     +  (Iopt_sgf
     +  ,y,z,x
     +  ,y0,z0,x0
     +  ,wall
     +  ,gyy,gyz,gyx
     +  ,gzy,gzz,gzx
     +  ,gxy,gxz,gxx
     +  ,presy,presz,presx
     +  ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +  ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +  ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +  )

      Txyy = Tyyx
      Txyz = Tzyx
      Tzyy = Tyyz

      Txzy = Tyzx
      Txzz = Tzzx
      Tzzy = Tyzz

      Txxy = Tyxx
      Txxz = Tzxx
      Tzxy = Tyxz

      www = wd(j)

      sumxx = sumxx - Txxy * www
      sumxy = sumxy - Txyy * www
      sumxz = sumxz - Txzy * www

      sumyx = sumyx - Tyxy * www
      sumyy = sumyy - Tyyy * www
      sumyz = sumyz - Tyzy * www

      sumzx = sumzx - Tzxy * www
      sumzy = sumzy - Tzyy * www
      sumzz = sumzz - Tzzy * www

      End Do

      sumxx = sumxx * rbases
      sumxy = sumxy * rbases
      sumxz = sumxz * rbases

      sumyx = sumyx * rbases
      sumyy = sumyy * rbases
      sumyz = sumyz * rbases

      sumzx = sumzx * rbases
      sumzy = sumzy * rbases
      sumzz = sumzz * rbases

      us = us - u0*sumxx - v0*sumyx  - w0*sumzx 
      vs = vs - u0*sumxy - v0*sumyy  - w0*sumzy 
      ws = ws - u0*sumxz - v0*sumyz  - w0*sumzz 

c     write (6,100)
c     write (6,100) i,sumxx,sumxy,sumxz
c     write (6,100) i,sumyx,sumyy,sumyz
c     write (6,100) i,sumzx,sumzy,sumzz
c     write (6,100)

c----------------  !   End of integration over the base

      us = us/pi4
      vs = vs/pi4
      ws = ws/pi4

      dlp(i,1) = us-u0
      dlp(i,2) = vs-v0
      dlp(i,3) = ws-w0

c     write (6,100) i,us,vs,ws

  1   Continue

c-----
c done
c-----

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      Return
      End
