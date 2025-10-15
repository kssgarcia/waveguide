      subroutine sgrad_3d_2p
     +
     +   (npts,nelm
     +   ,f
     +   ,gradfa
     +   )

c=======================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c=======================================

c----------------------------------
c computes the surface gradient
c of the scalar function f
c over a doubly-periodic surface
c
c SYMBOLS:
c --------
c
c components of grad(f);
c
c dfx:    df/dx
c dfy:    df/dy
c dfz:    df/dz

c dfxa:   df/dx averaged over host triangles
c dfya:   df/dy averaged over host triangles
c dfza:   df/dz averaged over host triangles
c
c  Iedge: interior/edge/corner node index:
c
c  Iedge(i,1) = 0 if global node i is an interior node
c
c  Iedge(i,1) = 1 if global node i is an edge node
c  Iedge(i,2): global index of image node
c
c  Iedge(i,1) = 3 if global node i is a corner node
c  Iedge(i,2): global index of first image node
c  Iedge(i,3): global index of first image node
c  Iedge(i,4): global index of first image node
c----------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension      p(1090,3)
      Dimension     ne(1090,7)
      Dimension    vna(1090,3)
      Dimension      f(1090)
      Dimension gradfa(1090,3),gradfa_sv(1090,3)
      Dimension itally(1090)

      Dimension Iedge(1090,4)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xxi(6),eet(6)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/edgepoints/Iedge

      common/albega/alpha,beta,gamma

      common/geo2/vna

c-----------
c initialize
c-----------

      Do i=1,npts
        gradfa(i,1) = 0.0
        gradfa(i,2) = 0.0
        gradfa(i,3) = 0.0
        itally(i) = 0
      End Do

c-------------------------------
c compute the gradient of f
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

      Do i=1,6      ! run over local nodes

        m = n(k,i)   ! global index of local node i on element k

        xi  = xxi(i)
        eta = eet(i)

        call sgrad_3d_2p_interp
     +
     +       (p(i1,1),p(i1,2),p(i1,3)
     +       ,p(i2,1),p(i2,2),p(i2,3)
     +       ,p(i3,1),p(i3,2),p(i3,3)
     +       ,p(i4,1),p(i4,2),p(i4,3)
     +       ,p(i5,1),p(i5,2),p(i5,3)
     +       ,p(i6,1),p(i6,2),p(i6,3)
     +       ,f(i1),f(i2),f(i3)
     +       ,f(i4),f(i5),f(i6)
     +       ,al,be,ga
     +       ,xi,eta
     +       ,vna(m,1),vna(m,2),vna(m,3)
     +       ,dfx,dfy,dfz
     +       )

        gradfa(m,1) = gradfa(m,1) + dfx
        gradfa(m,2) = gradfa(m,2) + dfy
        gradfa(m,3) = gradfa(m,3) + dfz

        itally(m) = itally(m)+1

      End Do

      End Do               ! loop over nodes

c----------------------------------
c average the tangential derivative
c over all elements
c-----------------------------------

      Do i = 1,npts
        par = float(itally(i))
        gradfa(i,1) = gradfa(i,1)/par
        gradfa(i,2) = gradfa(i,2)/par
        gradfa(i,3) = gradfa(i,3)/par
      End Do

c-----------------------
c account for the images
c-----------------------

c---
c save the points that have images
c---

      Do 76 i=1,npts

       if(Iedge(i,1).eq.0) Go to 76  ! no images

       gradfa_sv(i,1) = gradfa(i,1)
       gradfa_sv(i,2) = gradfa(i,2)
       gradfa_sv(i,3) = gradfa(i,3)

   76 Continue

c---
c average over the images
c---

      Do 77 i=1,npts

       noim = Iedge(i,1)   ! number of images

       If(noim.eq.0) Go to 77  ! no images

        Do j=2,noim+1     ! Iedge(i,1) images
         image = Iedge(i,j)
         gradfa(i,1) = gradfa_sv(i,1) + gradfa_sv(image,1)
         gradfa(i,2) = gradfa_sv(i,2) + gradfa_sv(image,2)
         gradfa(i,3) = gradfa_sv(i,3) + gradfa_sv(image,3)
        End Do

        gradfa(i,1) = gradfa(i,1)/noim
        gradfa(i,2) = gradfa(i,2)/noim
        gradfa(i,3) = gradfa(i,3)/noim

  77  Continue

c---
c done
c---

 100  Format (1x,i3,10(f12.8))
 101  Format (f12.8)

      return
      end
