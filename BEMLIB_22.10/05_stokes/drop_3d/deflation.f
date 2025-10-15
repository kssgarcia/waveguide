      subroutine  deflation 
     +
     +   (Idfl
     +   ,npts
     +   ,nelm
     +   ,mint
     +   ,vsk
     +   ,u
     +   ,CF1
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

c---------------------------------------------------
c Compute the projection coefficients for deflation
c
c NOTATION:
c
c vsk = (1-lamda)/(1+lamda)
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)
      Dimension   u(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/trq/xiq,etq,wq
      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c--------
c prepare
c--------

      SRFAR = 0.0D0  ! to compute surface area
      PRJ1  = 0.0D0  ! to compute u-dot-n

c---
c projections
c---

      Do k=1,nelm

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be =  beta(k)
      ga = gamma(k)

      SUM0 = 0.0D0   ! for surface area
      SUM1 = 0.0D0   ! for u-dot-n

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
     +   ,x,y,z                 ! output
     +   ,vnx,vny,vnz           ! output
     +   ,ux,uy,uz              ! output
     +   ,hs                    ! output
     +   )

       cf = 0.5D0*hs*wq(i)

       SUM0 = SUM0 + cf
       SUM1 = SUM1 + (vnx*ux + vny*uy + vnz*uz)*cf

       End Do

       SRFAR = SRFAR + SUM0
       PRJ1  = PRJ1  + SUM1

       End Do
 
c----------------------
c Deflation coefficient
c----------------------

       CF1 = vsk*PRJ1/SRFAR

c-----
c done
c-----

       return
       end
