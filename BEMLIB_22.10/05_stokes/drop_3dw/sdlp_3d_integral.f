      subroutine sdlp_3d_integral 
     +
     +  (x0,y0,z0
     +  ,u0,v0,w0
     +  ,k
     +  ,mint
     +  ,u
     +  ,uxel,uyel,uzel
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c---------------------------------------
c integrate the dlp over the kth element
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)
      Dimension   u(1026,3)

      Dimension     n(512,6),nbe (512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/var/shrt,wall

      common/trq/xiq,etq,wq

c------
c flags
c------

      Iopt_sgf = 2    ! for the Green's function

c----------------
c  prepare to run
c----------------

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be =  beta(k)
      ga = gamma(k)

      uxel = 0.0D0
      uyel = 0.0D0
      uzel = 0.0D0

c---
c triangle quadrature
c---

      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        call sdlp_3d_interp 
     +
     +   (p(i1,1),p(i1,2),p(i1,3)
     +   ,p(i2,1),p(i2,2),p(i2,3)
     +   ,p(i3,1),p(i3,2),p(i3,3)
     +   ,p(i4,1),p(i4,2),p(i4,3)
     +   ,p(i5,1),p(i5,2),p(i5,3)
     +   ,p(i6,1),p(i6,2),p(i6,3)
     +
     +   ,vna(i1,1),vna(i1,2),vna(i1,3)
     +   ,vna(i2,1),vna(i2,2),vna(i2,3)
     +   ,vna(i3,1),vna(i3,2),vna(i3,3)
     +   ,vna(i4,1),vna(i4,2),vna(i4,3)
     +   ,vna(i5,1),vna(i5,2),vna(i5,3)
     +   ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +   ,u(i1,1),u(i1,2),u(i1,3)
     +   ,u(i2,1),u(i2,2),u(i2,3)
     +   ,u(i3,1),u(i3,2),u(i3,3)
     +   ,u(i4,1),u(i4,2),u(i4,3)
     +   ,u(i5,1),u(i5,2),u(i5,3)
     +   ,u(i6,1),u(i6,2),u(i6,3)
     +
     +   ,al,be,ga
     +   ,xi,eta
     +
     +   ,x,y,z
     +   ,vnx,vny,vnz
     +   ,ux,uy,uz
     +   ,hs
     +   )

c---
c compute the Green's function
c---

      call sgf_3d_w 
     +
     +  (Iopt_sgf
     +  ,y,z,x
     +  ,y0,z0,x0
     +  ,wall
     +  ,gyy,gyz,gyx
     +  ,gzy,gzz,gzx
     +  ,gxy,gxz,gxx
     +  ,py,pz,px
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

      ux = ux-u0
      uy = uy-v0
      uz = uz-w0

      cf = 0.5D0*hs*wq(i)

c     ux = 0.0    ! for debugging
c     uy = 1.0
c     uz = 0.0

      uxel = uxel +((ux*Txxx + uy*Tyxx + uz*Tzxx)*vnx
     +             +(ux*Txxy + uy*Tyxy + uz*Tzxy)*vny
     +             +(ux*Txxz + uy*Tyxz + uz*Tzxz)*vnz)*cf

      uyel = uyel +((ux*Txyx + uy*Tyyx + uz*Tzyx)*vnx
     +             +(ux*Txyy + uy*Tyyy + uz*Tzyy)*vny
     +             +(ux*Txyz + uy*Tyyz + uz*Tzyz)*vnz)*cf

      uzel = uzel +((ux*Txzx + uy*Tyzx + uz*Tzzx)*vnx
     +             +(ux*Txzy + uy*Tyzy + uz*Tzzy)*vny
     +             +(ux*Txzz + uy*Tyzz + uz*Tzzz)*vnz)*cf

      End Do

c-----
c done
c-----

      return
      end
