      subroutine sgrad_3d
     +
     +    (npts,nelm
     +    ,f
     +    ,gradfa
     +    )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------
c Compute the surface gradient
c of a scalar function, f,
c defined at the nodes of a triangular grid
c over a three-dimensional surface
c
c SYMBOLS:
c --------
c
c dfx:  x component of grad(f) at the a node
c dfy:  y component of grad(f) at the a node
c dfz:  z component of grad(f) at the a node
c
c gradfa(i,j):  components of grad(f) at the ith node
c               averaged over host triangles: j=1,2,3
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension      p(1026,3)
      Dimension     ne(1026,7)
      Dimension    vna(1026,3)
      Dimension      f(1026)
      Dimension gradfa(1026,3)
      Dimension itally(1026)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xxi(6),eet(6)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo2/vna

c-----------
c initialize
c-----------

      Do i=1,npts

        gradfa(i,1) = 0.0D0
        gradfa(i,2) = 0.0D0
        gradfa(i,3) = 0.0D0

        Itally(i) = 0

      End Do

c-------------------------------
c compute the surface gradient of f
c at the nodes of each triangle
c-------------------------------

      Do k=1,nelm     ! run over elements

      i1 = n(k,1)
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

      Do i=1,6       ! run over local nodes

        m = n(k,i)   ! global index of local node i on element k

        xi  = xxi(i)
        eta = eet(i)

        call sgrad_3d_interp
     +
     +     (p(i1,1),p(i1,2),p(i1,3)
     +     ,p(i2,1),p(i2,2),p(i2,3)
     +     ,p(i3,1),p(i3,2),p(i3,3)
     +     ,p(i4,1),p(i4,2),p(i4,3)
     +     ,p(i5,1),p(i5,2),p(i5,3)
     +     ,p(i6,1),p(i6,2),p(i6,3)
     +
     +     ,f(i1),f(i2),f(i3)
     +     ,f(i4),f(i5),f(i6)
     +
     +     ,al,be,ga
     +     ,xi,eta
     +     ,vna(m,1),vna(m,2),vna(m,3)
     +     ,dfx,dfy,dfz
     +     )

        gradfa(m,1) = gradfa(m,1) + dfx
        gradfa(m,2) = gradfa(m,2) + dfy
        gradfa(m,3) = gradfa(m,3) + dfz

        Itally(m) = Itally(m)+1

      End Do

      End Do ! loop over nodes

c------------------------------------------
c average the tangential derivative
c over all elements
c------------------------------------------

      Do i=1,npts

        par = float(Itally(i))

        gradfa(i,1) = gradfa(i,1)/par
        gradfa(i,2) = gradfa(i,2)/par
        gradfa(i,3) = gradfa(i,3)/par

      End Do

c-----
c done
c-----

  99  Continue

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
