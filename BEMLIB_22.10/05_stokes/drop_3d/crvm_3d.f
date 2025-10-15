      subroutine crvm_3d
     +
     +   (nelm
     +   ,npts
     +   )

c=====================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=====================================

c----------------------------------------
c Compute the mean curvature at the 
c local element nodes,
c and then average to obtain the
c curvature at the global nodes
c
c  SYMBOLS:
c  --------
c
c  Itally(i): number of elements sharing a node
c             used for averaging
c
c  xxi, eet: local triangle coordinates
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension      p(1026,3)
      Dimension     ne(1026,7)
      Dimension    vna(1026,3)
      Dimension   crvm(1026)

      Dimension Itally(1026)    ! internal

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xxi(6),eet(6)   ! internal

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo2/vna
      common/geo3/crvm

c-----------
c initialize
c-----------

      Do i=1,npts
         crvm(i) = 0.0D0
       itally(i) = 0      ! for averaging at a node
      End Do

c----------------------
c loop over local nodes
c----------------------

      Do k=1,nelm    ! run over elements

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
c of the element nodes
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

       Do i=1,6  ! run over triangle nodes

        xi  = xxi(i)
        eta = eet(i)

        call crvm_3d_interp
     +
     +     (p(i1,1),p(i1,2),p(i1,3)
     +     ,p(i2,1),p(i2,2),p(i2,3)
     +     ,p(i3,1),p(i3,2),p(i3,3)
     +     ,p(i4,1),p(i4,2),p(i4,3)
     +     ,p(i5,1),p(i5,2),p(i5,3)
     +     ,p(i6,1),p(i6,2),p(i6,3)
     +
     +     ,vna(i1,1),vna(i1,2),vna(i1,3)
     +     ,vna(i2,1),vna(i2,2),vna(i2,3)
     +     ,vna(i3,1),vna(i3,2),vna(i3,3)
     +     ,vna(i4,1),vna(i4,2),vna(i4,3)
     +     ,vna(i5,1),vna(i5,2),vna(i5,3)
     +     ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +     ,al,be,ga
     +     ,xi,eta
     +     ,crvm_node
     +     )

        m = n(k,i)   ! global index of local node i
                     ! on element k

        crvm(m) = crvm(m) + crvm_node

        Itally(m) = Itally(m)+1

      End Do

      End Do       ! end of loop over elements

c---------------------------
c average the mean curvature
c at the nodes over elements
c---------------------------

      Do i=1,npts

        crvm(i) = crvm(i)/float(Itally(i))

      End Do

c-----
c done
c-----

  100 Format(1x,i3,10(1x,f10.6))
  101 Format(1x,i3,1x,i3,10(1x,f10.6))

      return
      end
