      subroutine vs_3d_circ_vort
     +
     +  (npts,nelm
     +  ,f
     +  ,vorta
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c----------------------------------
c Compute the vorticity
c at the nodes of a triangular grid
c by differentiating the circulation
c in the tangential plane
c
c SYMBOLS:
c -------
c
c vorta: vorticity at the nodes

c vorta: vorticity at the nodes
c        averaged over the neighboring elements
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)

      Dimension f(1026),vort(1026,3),vorta(1026,3)

      Dimension itally(1026)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xxi(6),eet(6)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma
      common/geo2/vna

c---
c initialize
c---

      Do i=1,npts

        vorta(i,1) = 0.0D0
        vorta(i,2) = 0.0D0
        vorta(i,3) = 0.0D0

        itally(i) = 0

      End Do

c-------------------------------
c compute the vorticity at the nodes
c of each triangle
c-------------------------------

      Do k=1,nelm     ! run over elements

       i1 = n(k,1)
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

       Do i=1,6    ! run over element nodes

        xi  = xxi(i)
        eta = eet(i)

        m = n(k,i)   ! global index of local node i on element k

        call vs_3d_circ_vort_interp
     +
     +      (p(i1,1),p(i1,2),p(i1,3)
     +      ,p(i2,1),p(i2,2),p(i2,3)
     +      ,p(i3,1),p(i3,2),p(i3,3)
     +      ,p(i4,1),p(i4,2),p(i4,3)
     +      ,p(i5,1),p(i5,2),p(i5,3)
     +      ,p(i6,1),p(i6,2),p(i6,3)
     +
     +      ,f(i1),f(i2),f(i3)
     +      ,f(i4),f(i5),f(i6)
     +
     +      ,al,be,ga
     +      ,xi,eta
     +
     +      ,vort(i,1),vort(i,2),vort(i,3)
     +      )

        vorta(m,1) = vorta(m,1)+vort(i,1)
        vorta(m,2) = vorta(m,2)+vort(i,2)
        vorta(m,3) = vorta(m,3)+vort(i,3)

        itally(m) = itally(m)+1

       End Do

      End Do  ! loop over nodes

c----------------------
c average the vorticity
c over all elements
c----------------------

      Do i=1,npts

        par = float(itally(i))

        vorta(i,1) = vorta(i,1)/par
        vorta(i,2) = vorta(i,2)/par
        vorta(i,3) = vorta(i,3)/par

      End Do

c----------------------------------
c project onto the tangential plane
c----------------------------------

      Do i=1,npts

       prj = vna(i,1)*vorta(i,1)
     +     + vna(i,2)*vorta(i,2)
     +     + vna(i,3)*vorta(i,3)

       vorta(i,1) = vorta(i,1)-prj*vna(i,1)
       vorta(i,2) = vorta(i,2)-prj*vna(i,2)
       vorta(i,3) = vorta(i,3)-prj*vna(i,3)

      End Do

c-----
c Done
c-----

  99  Continue

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
