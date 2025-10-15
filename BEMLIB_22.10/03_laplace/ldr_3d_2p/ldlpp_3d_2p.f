      subroutine ldlpp_3d_2p
     +
     +    (npts,nelm
     +    ,mint
     +    ,q
     +    ,x0,y0,z0
     +    ,ptl
     +    )

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement
c==========================================

c----------------------------------
c Compute the double-layer potential
c of a scalar function q 
c at the point (x0, y0, z0)
c
c SYMBOLS:
c -------
c
c ptl: double-layer potential at (x0, y0, z0)
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1090,3)
      Dimension   ne(1090,7)
      Dimension  vna(1090,3)

      Dimension  q(1090)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

      Parameter (tol=0.00000001)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

c---
c initialize
c---

      ptl = 0.0D0

c----------------------
c Compile the integrals
c over the triangles
c---------------------

       Do k=1,nelm     ! run over elements

c      write (6,*) " Integrating over element",k

       i1 = n(k,1)
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

c---
c apply the quadrature
c---

        call ldlpp_3d_2p_integral
     +
     +      (x0,y0,z0
     +      ,i,k
     +      ,mint
     +      ,q
     +      ,pptl
     +      ,arelm
     +      )

        ptl = ptl+pptl

        srf_area = srf_area+arelm

      End Do

c---
c done
c---

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
