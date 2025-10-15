      program trgl6_sqr_dr

c==========================================
c FDLIB  BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c==========================================

c-----------------------------------------------
c  Driver for the triangulation of a square
c  in the xy plane confined within -0.5<x,y<0.5
c  based on the successive subdivision of
c  and eight-element pattern
c
c  This program generates a structured grid
c  of six-node triangular elements
c
c  The z coordinate of the nodes can be adjusted
c  to give a doubly-periodic surface
c  with a desired geometry
c
c  SYMBOLS:
c  -------
c
c  p(i,j): coordinates of point i (j=1,2,3)
c
c  n(k,i): connectivity table: 
c          global labels of points for element k
c          i = 1, ..., 6
c
c  nbe(k,j): the three neighboring elements of element k (j=1,2,3)
c
c  ne(k,j):  ne(k,1) is the number of elements adjacent to point k
c            ne(k,2), ... are the elements numbers, j = 2, ..., 7
c            for this triangulation, up to six
c
c  npts:  total number of points
c  nelm:  total number of elements
c
c Iedge: interior-edge-corner node index:
c
c     Iedge(i,1) = 0 if global node i is an interior node
c
c     Iedge(i,1) = 1 if global node i is an edge node
c     Iedge(i,2): global index of image node
c
c     Iedge(i,1) = 3 if global node i is a corner node
c     Iedge(i,2): global index of first image node
c     Iedge(i,3): global index of first image node
c     Iedge(i,4): global index of first image node
c
c  alpha, beta, gamma: triangle mapping coefficients
c
c  vna: averaged normal vector at the nodes
c
c  arel: surface area of the elements
c
c  crvmel: mean curvature over an element
c         computed by contour integration
c
c  crvm: mean curvature at the nodes
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension     p(1090,3)
      Dimension    ne(1090,9)
      Dimension   vna(1090,3)
      Dimension  crvm(1090)
      Dimension Iedge(1090,4)
      Dimension   fnc(1090)
      Dimension color(1090)

      Dimension       n(512,6),nbe(512,3)
      Dimension   alpha(512),beta(512),gamma(512)
      Dimension    arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension  crvmel(512)

      Dimension xiq(20),etq(20),wq(20)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/edgepoints/Iedge

      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna
      common/geo3/crvm
      common/geo6/crvmel

      common/trq/xiq,etq,wq

c----------
c constants
c----------

      pi  = 3.14159265358 D0
      pi2 = 2.0D0*pi

      null = 0
      nfour = 4
      nseven = 7

c-----------
c input data
c-----------

  91  Continue

      write (6,*)
      write (6,*) " Will discretize a square patch"
      write (6,*) " Select the level of triangulation"
      write (6,*)
      write (6,*) " Choose from 0,1,2,3"
      write (6,*)
      write (6,*) " Enter 99 to quit"
      write (6,*) " ----------------"
      read  (5,*) ndiv 

      if(ndiv.eq.99) Go to 99

      if(ndiv.gt.3)  then
        write (6,*) 'out of range; please try again.'
        Go to 91
      end if

      RLx = 1.0D0
      RLy = 1.0D0

      write (6,*)
      write (6,*) " The z-elevation will be: "
      write (6,*) " z = ampx*cos(2 pi x/Lx) + ampy*cos(2 pi x/Ly)"
      write (6,*) " Please enter ampx and ampy"
      write (6,*) " --------------------------"
      read  (5,*) ampx,ampy

c----------------------
c triangulate the patch
c----------------------

      call trgl6_sqr (ndiv,npts,nelm)
      
      write (6,*)
      write (6,*) "Number of points:   ",npts
      write (6,*) "Number of elements: ",nelm
      write (6,*)

c---
c set the wave numbers
c---

      wnx = pi2/RLx
      wny = pi2/RLy

c---
c set the z coordinate
c of global nodes
c---


      Do i=1,npts
        p(i,3) = ampx*Dcos(wnx*p(i,1))
     +          +ampy*Dcos(wny*p(i,2))
      End Do

c----
c print the edge nodes
c       and their images
c---

c     write (6,*)
c     write (6,*) " Edge nodes and images"
c     write (6,*)

c     Do i=1,npts
c      If(Iedge(i,1).ne.0) then
c       write (6,139) i,(Iedge(i,j),j=1,1+Iedge(i,1))
c      End If
c     End Do

      strain = 0.0;

c-----------------
c open output file
c-----------------

c-----------------------------

 98   Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 to see the whole surface"
      write (6,*) " 2 to see a particular element"
      write (6,*) " 3 to compute the surface area and volume"
      write (6,*) " 4 to compute the surface area, volume,"
      write (6,*) "   normal vector, and mean curvature"
      write (6,*) " 5 to compute the surface integral"
      write (6,*) "    of a function"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c-----------------------
      if(menu.eq.1) then
c-----------------------

      open (1,file="trgl_sqr.net")

c print a matlab file

      Do i=1,npts
         color(i)=p(i,3)
      End Do

      Index = 2  ! 3-node triangles arising by subdivision
      Index = 1  ! 6-node triangles

      if(Index.eq.1) then
        write (1,*) nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      elseif(Index.eq.2) then
        write (1,*) nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end if

      write (1,*) strain

      Do k=1,nelm
        write (1,*) 
        call printel(k,Index,color)
      End Do

      write (1,100) null
      close (1)

      Go to 98

      end if

c----------------------------
      if(menu.eq.2) then
c----------------------------

      open (1,file="trgl_sqr.net")

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " the number of a triangle"
      write (6,*) " to see it and its neighbors"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) k

      if(k.eq.0) Go to 99

      Do i=1,npts
        color(i)=p(i,3)
      End Do

      Index = 2  ! 3-node triangles arising by subdivision
      Index = 1  ! 6-node triangles

      if(Index.eq.1) then
        write (1,*) nseven
        write (1,*) 7*nelm
        write (1,*) 4
      elseif(Index.eq.2) then
        write (1,*) nfour
        write (1,*) 16*4
        write (1,*) 4*4
      end if

      write (1,*) strain

c---
c print the element and its three neighbors
c---

      call printel (k,Index,color)
      call printel (nbe(k,1),Index,color)
      call printel (nbe(k,2),Index,color)
      call printel (nbe(k,3),Index,color)

      write (1,100) null
      close (1)

      Go to 98

      End If

c---------------------------------------
c Will perform various computations
c based on the parametric representation
c---------------------------------------

c---
c Select triangle quadrature
c for surface integration
c---

       write (6,*)
       write (6,*) " Will use the m-point gauss-triangle rule"
       write (6,*) "       for integration over each triangle"
       write (6,*) " Please enter m"
       write (6,*) " choose from 1,3,4,6,7,9,12,13"
       write (6,*) " Enter 0 to quit"
       write (6,*) "----------------"
       read  (5,*) mint

       If(mint.eq.0) Go to 99

       call gauss_trgl (mint,xiq,etq,wq)

c---
c Compute the coefficients:
c alpha, beta, gamma,
c for the xi-eta mapping
c over each triangle
c---

      Do k=1,nelm

       i1 = n(k,1)
       i2 = n(k,2)
       i3 = n(k,3)
       i4 = n(k,4)
       i5 = n(k,5)
       i6 = n(k,6)

       call abc
     +
     +   (p(i1,1),p(i1,2),p(i1,3)
     +   ,p(i2,1),p(i2,2),p(i2,3)
     +   ,p(i3,1),p(i3,2),p(i3,3)
     +   ,p(i4,1),p(i4,2),p(i4,3)
     +   ,p(i5,1),p(i5,2),p(i5,3)
     +   ,p(i6,1),p(i6,2),p(i6,3)
     +   ,alpha(k),beta(k),gamma(k)
     +   )

       End Do

c-----------------------
      If(menu.eq.3) then
c-----------------------

c------------------------------------
c Compute the area and the volume
c------------------------------------

c---
c initialize
c---

      area = 0.0
      vlm  = 0.0

c---
c loop over the elements
c---

      Do k=1,nelm

        i1 = n(k,1)
        i2 = n(k,2)
        i3 = n(k,3)
        i4 = n(k,4)
        i5 = n(k,5)
        i6 = n(k,6)

c---
c loop over the integration points
c---

        Do i=1,mint

          xi  = xiq(i)
          eta = etq(i)

c---
c interpolate for the normal vector
c and surface metric
c---

          call interp_vn_hs
     +
     +      (p(i1,1),p(i1,2),p(i1,3)
     +      ,p(i2,1),p(i2,2),p(i2,3)
     +      ,p(i3,1),p(i3,2),p(i3,3)
     +      ,p(i4,1),p(i4,2),p(i4,3)
     +      ,p(i5,1),p(i5,2),p(i5,3)
     +      ,p(i6,1),p(i6,2),p(i6,3)
     +      ,alpha(k),beta(k),gamma(k)
     +      ,xi,eta
     +      ,x,y,z
     +      ,vnx,vny,vnz
     +      ,hs
     +      )

          area = area +       hs*wq(i)
          vlm  = vlm  + z*vnz*hs*wq(i)

        End Do

      End Do

      area = 0.5*area     ! factor 0.5 from the quadarature
      vlm  = 0.5*vlm      ! factor 0.5 from the quadarature

      write (6,*)
      write (6,110) area
      write (6,111) vlm
      write (6,*)

c----------------------------
      elseif(menu.eq.4) then
c----------------------------

      open (1,file="trgl_sqr.net")

c------
c Compute:
c
c (a) the surface area of the individual elements.
c (b) the x, y, and z surface moments over each element.
c (c) the total surface area and volume.
c (d) the mean curvature of each element.
c (e) the average value of the normal vector at each node;
c     this is done by computing the normal vector at the 6 nodes
c     of each triangle, and then averaging the contributions.
c------

      call elm_geom
     +
     +   (nelm,npts,mint
     +   ,xmom,ymom,zmom
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c---
c display
c---

      write (6,*)
      write (6,110) area
      write (6,111) vlm
      write (6,*)

      write (6,*)
      write (6,*) " Element surface area and mean curvature"
      write (6,*)

      Do i=1,nelm
       write (6,100) i,arel(i),crvmel(i)
      End Do

c------
c Compute the averaged mean curvature at the nodes
c
c Note that subroutine: crvm_3d_2p
c should be called after: elm_geom
c------

      call crvm_3d_2p (nelm,npts)

      write (6,*)
      write (6,*) " Node position, normal vector, mean curvature"
      write (6,*)

      Do i=1,npts
       write (6,102) i,p(i,1),p(i,2),p(i,3)
     +                ,vna(i,1),vna(i,2),vna(i,3)
     +                ,crvm(i)
      End Do

      Do i=1,npts
         color(i)=crvm(i)
      End Do

      Index = 2

      if(Index.eq.1) then
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      elseif(Index.eq.2) then
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end If

      write (1,*) strain

      Do k=1,nelm
        call printel(k,Index,color)
      End Do

      write (1,100) null
      close (1)

c----------------------------
      Else If(menu.eq.5) then
c----------------------------

c----------------------------
c Compute the surface integral
c of a function
c----------------------------

 93   Continue

      write (6,*)
      write (6,*) " MENU OF FUNCTIONS"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 for f = 1.0"
      write (6,*) "       2 for f = z**2"
      write (6,*) "       3 for f = z**3"
      write (6,*) "---------------------"
      read  (5,*) my_choice

      If(my_choice.eq.0) Go to 98

c---
c evaluate the integrand at the nodes
c---

      Do i=1,npts
        If(my_choice.eq.1) fnc(i) = 1.0
        If(my_choice.eq.2) fnc(i) = p(i,3)**2
        If(my_choice.eq.3) fnc(i) = p(i,3)**3
      End Do

c---
c perform the integration
c---

      call srf_int_3d_2p
     +
     +  (nelm,npts,mint
     +  ,fnc
     +  ,olok
     +  )

      write (6,*)
      write (6,456) olok

      Go to 93

c-----------
      End If
c-----------

      Go to 98    ! return to the menu

c-------
c Ending
c-------

  99  Continue

c-----
c Done
c-----

  100 Format (1x,i4,10(1x,f12.5))
  101 Format (10(1x,f12.5))
  102 Format (1x,i4,10(1x,f8.5))
  103 Format (10(1x,f8.5))

  110 Format (" Surface Area :",F15.10)
  111 Format (" Volume       :",F15.10)

  456 Format ("Integral = ",f15.10)

  139 Format (i4," : ",i1,";",9(1x,i4))

      Stop
      End
