      subroutine srf_int_3d
     +
     +  (nelm,npts
     +  ,mint
     +  ,fnc
     +  ,olok
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c---------------------------------------
c Compute:
c
c the surface integral of a function fnc
c using a quadrature
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension fnc(1026)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension  xxi(6), eet(6)

      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/trq/xiq,etq,wq

c-----------
c initialize
c-----------

      olok = 0.0

c------------------------
c perform the integration
c------------------------

      Do k=1,nelm    ! loop over elements

       i1 = n(k,1)     ! global node labels
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       al = alpha(k)
       be = beta (k)
       ga = gamma(k)

c---
c triangle coordinates
c of the nodes
c---

       xxi(1) = 0.0D0
       eet(1) = 0.0D0

       xxi(2) = 1.0D0
       eet(2) = 0.0D0

       xxi(3) = 0.0D0
       eet(3) = 1.0D0

       xxi(4) = al
       eet(4) = 0.0D0

       xxi(5) = ga
       eet(5) = 1.0D0-ga

       xxi(6) = 0.0D0
       eet(6) = be

c---
c quadrature
c---

       Do i=1,mint    ! loop over integration points

        xi  = xiq(i)
        eta = etq(i)

        call srf_int_3d_interp
     +
     +     (p(i1,1),p(i1,2),p(i1,3)
     +     ,p(i2,1),p(i2,2),p(i2,3)
     +     ,p(i3,1),p(i3,2),p(i3,3)
     +     ,p(i4,1),p(i4,2),p(i4,3)
     +     ,p(i5,1),p(i5,2),p(i5,3)
     +     ,p(i6,1),p(i6,2),p(i6,3)
     +     ,fnc(i1),fnc(i2),fnc(i3)
     +     ,fnc(i4),fnc(i5),fnc(i6)
     +
     +     ,al,be,ga
     +     ,xi,eta
     +     ,x,y,z
     +
     +     ,hs
     +
     +     ,f
     +     )

        olok = olok + f*hs*wq(i)

       End Do

      End Do ! End of loop over elements

      olok = 0.5D0*olok  ! factor 0.5 due to the quadrature

c-----
c done
c-----

      return
      end
