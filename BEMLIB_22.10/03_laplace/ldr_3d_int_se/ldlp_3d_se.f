      subroutine ldlp_3d_se
     +
     +   (npts
     +   ,nelm
     +   ,intm   ! integration method
     +   ,mint
     +   ,nsp    ! spectral connectivity
     +   ,psp    ! spectral global nodes
     +   ,q      ! spectral dipole
     +   ,qsint  !   surface integral
     +   ,dlp    ! spectral double layer
     +   )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c----------------------------------
c Compute the principal value of the
c double-layer potential
c of a scalar function q
c at the nodes of a triangular grid
c on a closed surface
c
c Also compute the surface integral of q
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c spectral:

      Dimension   nsp(512,100)   ! connectivity
      Dimension   psp(20000,3)
      Dimension   q(20000),dlp(20000)

      Parameter (tol=0.00000001)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

      common/trq/xiq,etq,wq

      common/spectral2/nesp,ngsp

c----------------------
c launch the quadrature
c----------------------

      Do i=1,ngsp   ! loop over spectral nodes

c      write (6,*) " Computing the dlp potential at point",i

       x0 = psp(i,1)
       y0 = psp(i,2)
       z0 = psp(i,3)

       q0 = q(i)

c      write (6,100) i,x0,y0,z0,q0

c-----------
c initialize
c-----------

       srf_area = 0.0D0
       ptl      = 0.0D0
       qsint    = 0.0D0

c----------------------
c compile the integrals
c over the triangles
c---------------------

       Do k=1,nelm     ! run over elements

c      write (6,*) " Integrating over element",k

c---
c apply the quadrature
c---

        call ldlp_3d_se_integral
     +
     +    (x0,y0,z0
     +    ,nsp
     +    ,psp
     +    ,k
     +    ,intm
     +    ,mint
     +    ,q
     +    ,q0
     +    ,pptl
     +    ,aarel
     +    ,qqsint
     +    )

c       write (6,*) pptl,aarel,qqsint

        ptl = ptl+pptl

        srf_area = srf_area + aarel
        qsint    = qsint + qqsint

      End Do

c--------------------
c      If(i.eq.7) then
c       write (6,*)
c       write (6,*) " total surface area computed in ldlp: ",srf_area
c       write (6,*)
c      End If
c--------------------

c---
c account for the principal value
c---

       dlp(i) = ptl - 0.50D0 * q0
c      write(6,*) dlp(i),qsint,q0

      End Do               ! loop over nodes

c-----
c done
c-----

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
