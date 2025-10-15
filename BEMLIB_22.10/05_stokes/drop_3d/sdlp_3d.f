      subroutine sdlp_3d 
     +
     +   (npts
     +   ,nelm
     +   ,mint
     +   ,Iflow
     +   ,u
     +   ,dlp
     +   )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------------------------
c Compute the principal value of the double-layer
c potential at the surface nodes
c
c Desingularization is done by using an
c integral identity
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)
      Dimension   u(1026,3)

      Dimension nvel(1026)
      Dimension  dlp(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/var/shrt,wall

      common/veloc1/nvelt,nvel

      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c---
c for triply-periodic flow
c---

      common/aaaa/a11,a12,a13,a21,a22,a23,a31,a32,a33
      common/bbbb/b11,b12,b13,b21,b22,b23,b31,b32,b33
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

c----------------------
c     Do i=1,npts     ! to debug
c       u(i,1) = 0.   ! to debug
c       u(i,2) = 0.   ! to debug
c       u(i,3) = 0.   ! to debug
c     End Do          ! to debug
c---------------------

c----------------
c loop over nodes
c----------------

      Do node=1,nvelt

      i = nvel(node)      ! global node label

      x0 = p(i,1)
      y0 = p(i,2)
      z0 = p(i,3)

      u0 = u(i,1)
      v0 = u(i,2)
      w0 = u(i,3)

c--------------------------------------
c     write (6,*) " Enter x0,y0,z0"    ! to debug
c     read  (5,*) x0,y0,z0             ! to debug
c     write (6,*) " Enter u0,v0,w0"    ! to debug
c     read  (5,*) u0,v0,w0             ! to debug
c--------------------------------------

      us = 0.0D0
      vs = 0.0D0
      ws = 0.0D0

      Do k=1,nelm

        call sdlp_3d_integral 
     +
     +   (x0,y0,z0
     +   ,u0,v0,w0
     +   ,k
     +   ,mint
     +   ,Iflow
     +   ,u
     +   ,uxel,uyel,uzel
     +   )

        us = us + uxel
        vs = vs + uyel
        ws = ws + uzel

      End Do

      us = us/pi4
      vs = vs/pi4
      ws = ws/pi4

      dlp(i,1) = us-u0
      dlp(i,2) = vs-v0
      dlp(i,3) = ws-w0

c     write (6,100) i,us,vs,ws

      End Do

c-----
c done
c-----

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
