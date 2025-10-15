      subroutine srf_int_3d_2p
     +
     +                     (nelm
     +                     ,npts
     +                     ,mint
     +                     ,fnc
     +                     ,olok
     +                     )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c------------------------------------------
c Compute:
c
c the surface integral of the function: fnc
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1090,3)
      Dimension  ne(1090,7)
      Dimension fnc(1090)

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

      Do 1 k=1,nelm    ! loop over elements

       i1 = n(k,1)     ! global node labels
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       al = alpha(k)
       be =  beta(k)
       ga = gamma(k)

c---
c triangle coordinates
c of the nodes
c---

       xxi(1) = 0.0
       eet(1) = 0.0

       xxi(2) = 1.0
       eet(2) = 0.0

       xxi(3) = 0.0
       eet(3) = 1.0

       xxi(4) = al
       eet(4) = 0.0

       xxi(5) = ga
       eet(5) = 1.0-ga

       xxi(6) = 0.0
       eet(6) = be

c---
c quadrature
c---

       Do i = 1,mint    ! loop over integration points

        xi  = xiq(i)
        eta = etq(i)

        call interp_srf_int_3d_2p
     +
     +                (p(i1,1),p(i1,2),p(i1,3)
     +                ,p(i2,1),p(i2,2),p(i2,3)
     +                ,p(i3,1),p(i3,2),p(i3,3)
     +                ,p(i4,1),p(i4,2),p(i4,3)
     +                ,p(i5,1),p(i5,2),p(i5,3)
     +                ,p(i6,1),p(i6,2),p(i6,3)
     +                ,fnc(i1),fnc(i2),fnc(i3)
     +                ,fnc(i4),fnc(i5),fnc(i6)
     +                ,al,be,ga
     +                ,xi,eta
     +                ,x,y,z
     +                ,hs
     +                ,f
     +                )

        olok = olok + f*hs*wq(i)

       End Do

  1   Continue         ! End of look pver elements


      olok = 0.5*olok  ! factor 0.5 due to the quadrature

c---
c Done
c---

      Return
      End

c=================================================

      subroutine interp_srf_int_3d_2p
     +
     +                    (x1,y1,z1
     +                    ,x2,y2,z2
     +                    ,x3,y3,z3
     +                    ,x4,y4,z4
     +                    ,x5,y5,z5
     +                    ,x6,y6,z6
     +                    ,f1,f2,f3,f4,f5,f6
     +                    ,al,be,ga
     +                    ,xi,eta
     +                    ,x,y,z
     +                    ,hs
     +                    ,f
     +                    )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c--------------------------------
c Interpolates over an element
c to compute following variables:
c
c   position vector
c   surface metric
c   surface function f
c---------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c prepare
c---

      alc = 1.0-al
      bec = 1.0-be
      gac = 1.0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c---
c  interpolation basis functions
c---

      ph2 = xi *(xi -al+eta*(al-ga)/gac)   /alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0)/ga)/bec
      ph4 = xi *(1.0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0-xi-eta)/bebec
      ph1 = 1.0-ph2-ph3-ph4-ph5-ph6

c---
c  interpolate position vector (x, y, z)
c  and function f
c---

      x = x1*ph1 + x2*ph2 + x3*ph3 
     +  + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3
     +  + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3
     +  + z4*ph4 + z5*ph5 + z6*ph6

      f = f1*ph1 + f2*ph2 + f3*ph3
     +  + f4*ph4 + f5*ph5 + f6*ph6

c---
c  xi derivatives of phi
c---

      dph2 =  (2.0*xi-al +eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0)/(ga*bec)
      dph4 =  (1.0-2.0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6


c---
c  compute dx/dxi from xi derivatives of phi
c---

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4
     +      + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4
     +      + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4
     +      + z5*dph5 + z6*dph6

c---
c  eta derivatives of phi
c---

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0*eta-be+xi*(be+ga-1.0)/ga)/bec
      pph4 =  -xi/alalc
      pph5 =   xi/gagac
      pph6 =  (1.0-xi-2.0*eta)/bebec
      pph1 = - pph2-pph3-pph4-pph5-pph6

c---
c  compute dx/deta from eta derivatives of phi
c---

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4
     +      + x5*pph5 + x6*pph6
      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4
     +      + y5*pph5 + y6*pph6
      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4
     +      + z5*pph5 + z6*pph6

c---
c normal (non-unit) vector:
c vn = (DxDxi)x(DxDeta)
c
c surface metric:  hs = norm(vn)
c---

      vnx = DyDxi*DzDet - DyDet*DzDxi
      vny = DzDxi*DxDet - DzDet*DxDxi
      vnz = DxDxi*DyDet - DxDet*DyDxi

      hs  = sqrt(vnx**2+vny**2+vnz**2)

c---
c Done
c---

  99  Continue

      Return
      End
