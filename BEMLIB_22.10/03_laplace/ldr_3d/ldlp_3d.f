      subroutine ldlp_3d
     +
     +   (npts,nelm
     +   ,mint
     +   ,q
     +   ,Iflow
     +   ,dlp
     +   )

c===========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c----------------------------------
c Compute the principal value of the
c double-layer potential
c of a scalar function q
c at the nodes of a triangular grid
c on a closed surface
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)

      Dimension  q(1026),dlp(1026)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

      Parameter (tol=0.00000001)

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

      Do i=1,npts   ! loop over nodes

c      write (6,*) " Computing the dlp at point",i

       x0 = p(i,1)
       y0 = p(i,2)
       z0 = p(i,3)

       q0 = q(i)

c-----------
c initialize
c-----------

       srf_area = 0.0D0
       ptl      = 0.0D0

c----------------------
c Compile the integrals
c over the triangles
c---------------------

       Do k=1,nelm     ! run over elements

c      write (6,*) " lpdl_3d: ntegrating over element",k

c---
c integration will be performed 
c only if q is nonzero
c---

       i1 = n(k,1)
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       test = Dabs(q(i1))+Dabs(q(i2))+Dabs(q(i3))
     +      + Dabs(q(i4))+Dabs(q(i5))+Dabs(q(i6))
     +      + Dabs(q0)

c---
c apply the quadrature
c---

        if(test.gt.tol) then

        call ldlp_3d_integral
     +
     +    (x0,y0,z0
     +    ,i,k
     +    ,mint
     +    ,Iflow
     +    ,q
     +    ,q0
     +    ,pptl
     +    ,arelm
     +    )

        ptl = ptl + pptl

        srf_area = srf_area+arelm

        end if

       End Do

c--------------------
c      if(i.eq.7) then
c       write (6,*)
c       write (6,*) " total surface area computed in ldlp: ",srf_area
c       write (6,*)
c      End If
c--------------------

c---
c account for the principal value
c---

       dlp(i) = ptl - 0.50D0 * q0

      End Do    ! loop over nodes

c-----
c done
c-----

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      Return
      End
