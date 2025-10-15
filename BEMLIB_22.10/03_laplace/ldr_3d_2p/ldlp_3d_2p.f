      subroutine ldlp_3d_2p
     +
     +  (npts,nelm
     +  ,mint
     +  ,q
     +  ,dlp
     +  )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c----------------------------------------

c----------------------------------
c Compute the double-layer potential
c of a scalar function q
c at the nodes of a triangular grid
c on a doubly-periodic surface
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1090,3)
      Dimension   ne(1090,7)
      Dimension  vna(1090,3)

      Dimension  q(1090),dlp(1090)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

      Parameter (tol=0.000000001)

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

      Do i = 1,npts

c      write (6,*) " Computing the dlp potential at point",i

       x0 = p(i,1)
       y0 = p(i,2)
       z0 = p(i,3)

       q0 = q(i)

c-----
c initialize
c-----

       srf_area = 0.0D0
       ptl      = 0.0D0

c----------------------
c Compile the integrals
c over the triangles
c---------------------

       Do 2 k=1,nelm     ! run over elements

c      write (6,*) " Integrating over element",k

c---
c integration will be performed only
c if q is nonzero
c---

       i1 = n(k,1)
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       test = abs(q(i1))+abs(q(i2))+abs(q(i3))
     +      + abs(q(i4))+abs(q(i5))+abs(q(i6))
     +      + abs(q0)

       if(test.le.tol) Go to 2

c---
c apply the quadrature
c---

        call ldlp_3d_2p_integral
     +
     +     (x0,y0,z0
     +     ,i,k
     +     ,mint
     +     ,q,q0
     +     ,pptl
     +     ,arelm
     +     )

        ptl = ptl+pptl

        srf_area = srf_area+arelm

  2   Continue

c--------------------
c      if(i.eq.7) then
c       write (6,*)
c       write (6,*) " total surface area computed in ldlp: ",srf_area
c       write (6,*)
c      end if
c--------------------

c---
c account for the principal value
c---

       dlp(i) = ptl - 0.5*q0

       End Do               ! loop over nodes

c---
c Done
c---

  99  Continue

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      Return
      End
