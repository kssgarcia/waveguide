      subroutine vtg_3d
     +
     + (nelm
     + ,npts
     + ,v
     + ,vtg
     + )

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------------
c Compute the tangential gradient of a vector
c at the nodes: vtg
c
c This is done by interpolating and then averaging 
c over elements sharing each node
c
c vtg is computed by solving three linear systems
c for three equations for the columns
c arising from the equations:
c
c dv/dxi = dx/dxi . vtg
c dv/det = dx/det . vtg
c 0      = n      . vtg
c
c where n is the normal vector
c
c
c  SYMBOLS:
c  --------
c
c  vtg(i,j,k): tangential gradient of vector: v
c              averaged at the ith node
c
c  itally(i): number of elements sharing a node
c             for averaging
c
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    p(1026,3)
      Dimension   ne(1026,7)
      Dimension  vna(1026,3)

      Dimension    v(1026,3)
      Dimension  vtg(1026,3,3)
      Dimension vtgn(3,3)
      Dimension Itally(1026)     ! internal

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xxi(6), eet(6)

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

       Do j=1,3
        Do k=1,3
         vtg(i,j,k) = 0.0D0
        End Do
       End Do
 
       itally(i) = 0 

      End Do

c------------------

      Do 1 k=1,nelm    ! loop over elements

       i1 = n(k,1)     ! global node labels
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       al = alpha(k)
       be = beta (k)
       ga = gamma(k)

c---------------------
c triangle coordinates
c at the nodes
c---------------------

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

c---------------------------------------
c compute the vector tangential gradient
c---------------------------------------

       Do i=1,6    ! loop over element-nodes

        m = n(k,i)   ! global index of local node i on element k

        xi  = xxi(i)
        eta = eet(i)

        call vtg_3d_interp
     +
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +
     +    ,vna(i1,1),vna(i1,2),vna(i1,3)
     +    ,vna(i2,1),vna(i2,2),vna(i2,3)
     +    ,vna(i3,1),vna(i3,2),vna(i3,3)
     +    ,vna(i4,1),vna(i4,2),vna(i4,3)
     +    ,vna(i5,1),vna(i5,2),vna(i5,3)
     +    ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +    ,v(i1,1),v(i1,2),v(i1,3)
     +    ,v(i2,1),v(i2,2),v(i2,3)
     +    ,v(i3,1),v(i3,2),v(i3,3)
     +    ,v(i4,1),v(i4,2),v(i4,3)
     +    ,v(i5,1),v(i5,2),v(i5,3)
     +    ,v(i6,1),v(i6,2),v(i6,3)
     +
     +    ,al,be,ga
     +    ,xi,eta
     +    ,vtgn
     +    )

        Do j=1,3
         Do l=1,3
           vtg(m,j,l) = vtg(m,j,l) + vtgn(j,l)
          End Do
        End Do

        itally(m) = itally(m)+1

      End Do

  1   Continue       ! end of loop over elements

c---------------------------
c average the curvature
c tensor at the nodes
c---------------------------

      Do i=1,npts

        par = float(itally(i))

        Do j=1,3
         Do l=1,3
          vtg(i,j,l) = vtg(i,j,l)/par
         End Do
        End Do

      End Do

c-----
c Done
c-----

      return
      end
