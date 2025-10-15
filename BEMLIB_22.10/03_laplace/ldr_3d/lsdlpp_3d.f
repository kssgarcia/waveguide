      subroutine lsdlpp_3d
     +
     +    (nelm
     +    ,mint
     +    ,f
     +    ,dfdn
     +    ,x0,y0,z0
     +    ,Iflow
     +    ,f0
     +    )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c---------------------------------------
c Evaluate a harmonic function at a point
c using the boundary integral representation
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)

      Dimension  f(1026),dfdn(1026)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

c----------------------
c launch the quadrature
c----------------------

c-----
c initialize
c-----

       srf_area = 0.0D0
       f0       = 0.0D0

c----------------------
c compile the integrals
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

        call lsdlpp_3d_integral
     +
     +   (x0,y0,z0
     +   ,k
     +   ,mint,Iflow
     +   ,f
     _   ,dfdn
     +   ,sdlp
     +   ,arelm
     +   )

        srf_area = srf_area+arelm

        f0 = f0 + sdlp

      End Do

c--------------------
c      write (6,*)
c      write (6,*) " total surface area computed in lsdlpp: ",srf_area
c      write (6,*)
c--------------------

c-----
c done
c-----

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
