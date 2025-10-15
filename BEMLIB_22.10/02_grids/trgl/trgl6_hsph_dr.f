      program trgl6_hsph_dr

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c-----------------------------------
c  Driver for triangulating a section
c  of a sphere resting on y=0
c  into 6-node triangles.
c
c  Triangulation is carried out
c  by subdividing half a regular octahedron
c  into descendant elements
c
c  This program generates an ordered grid
c  of six-node triangular elements
c  on the surface of a sphere.
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
c  nbe(k,j):  the three neighboring elements of element k (j=1,2,3)

c  ne(k,j):  ne(k,1) is the number of elements adjacent to point k
c            ne(k,2), ... are the elements numbers, j = 2, ..., 7
c            for this triangulation, up to six
c
c  npts:  total number of points
c  nelm:   total number of elements
c
c  alpha, beta, gamma: triangle mapping coefficients
c
c  vna: averaged normal vector at the nodes
c
c  arel: surface area of the elements
c
c  crvmel: mean curvature over an element
c          computed by contour integration
c
c  crvm: mean curvature at the nodes
c-----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension     p(1026,3)
      Dimension    ne(1026,7)
      Dimension   vna(1026,3)
      Dimension  crvm(1026)
      Dimension color(1026)

      Dimension       n(512,6), nbe(512,3)
      Dimension   alpha(512),  beta(512),gamma(512)
      Dimension    arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension  crvmel(512)

      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna
      common/geo3/crvm
      common/geo6/crvmel

      common/trq/xiq,etq,wq

c----------
c constants
c----------

      pi  = 3.1415926 5358 D0
      pi2 = 2.0D0*pi
      pih = 0.5D0*pi

      Null   = 0
      Nfour  = 4
      Nseven = 7
      oot  = 1.0D0/3.0D0

c------------
c preferences
c------------

      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to read input from file: trgl6_hsph.dat"
      write (6,*) " 2 to enter data"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ienrd

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (1,file="trgl6_hsph.dat")

       read (1,*) Ndiv
       read (1,*)
       read (1,*) cont_angle
       read (1,*) req

      close (1)

c------------------------
      else if(Ienrd.eq.2) then
c------------------------

   91 write (6,*)
      write (6,*) " Select the level of triangulation"
      write (6,*)
      write (6,*) " Enter 0, 1, 2, 3"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) ndiv 

      if(ndiv.eq.99) Go to 99

      if(ndiv.gt.3) then
        write (6,*) "Out of range: please try again"
        Go to 91
      end if

      write (6,*)
      write (6,*) "Enter the contact angle in multiples of pi"
      write (6,*)
      write (6,*) "0.5 for a hemisphere"
      write (6,*) "1.0 for a whole sphere"
      write (6,*) "99 to quit"
      write (6,*) "-----------"
      read (5,*) cont_angle

      if(cont_angle.eq.99) Go to 99

      write (6,*)
      write (6,*) "The volume enclosed by the section of the"
      write (6,*) "sphere and its base will be equal to:"
      write (6,*)
      write (6,*) " 2*pi/3 * a**3,"
      write (6,*)
      write (6,*) "where a is the equivalent radius"
      write (6,*)
      write (6,*) "Please enter the equivalent radius: a"
      write (6,*) "-------------------------------------"
      write (6,*)
      read (5,*) req

c-----------
      end if
c-----------

c---------------------------
c triangulate the hemisphere
c---------------------------

      call trgl6_hsph_octa (ndiv,npts,nelm)
      
      write (6,*)
      write (6,*) "Number of points:   ",npts
      write (6,*) "Number of elements: ",nelm

c----------------------------------------
c  Perform the following transformations:
c
c  1) Slide the hemisphere to 
c     to the specified contact angle.
c
c  2) Translate along the y axis so that
c     the contact line is at y = 0.
c
c  3) Scale so that the volume of the
c     spherical section is identical
c     to that of the unit hemishpere.
c----------------------------------------

      phi   = pi*cont_angle
      cs    = Dcos(phi)

      scale = 2.0D0/3.0D0 * 1.0D0/(1.0D0-cs-(1.0D0-cs**3)/3.0D0)
      scale = req*scale**oot

c---- over-ride
      scale = 1.0D0
c---- over-ride

      shift = Dsin(pih-phi)

      Do i=1,npts

       rr = Dsqrt(p(i,1)**2+p(i,3)**2)

       thet0 = Datan2(rr,p(i,2))
       thetz = Datan2(p(i,1),p(i,3))

       thet = thet0*phi/pih

       rrr = Dsin(thet)

       p(i,1) = rrr*sin(thetz)
       p(i,2) =     cos(thet)
       p(i,3) = rrr*cos(thetz)

       p(i,2) = p(i,2)-shift

       p(i,1) = p(i,1)*scale
       p(i,2) = p(i,2)*scale
       p(i,3) = p(i,3)*scale

      End Do

c---
c define the surface color
c---

      Do i=1,npts
       color(i) = p(i,2)
      End Do

c------------------
c open output files
c------------------

      open (1,file="trgl6_hsph.net")
      open (8,file="trgl6_hsph.stl")
      write (8,*) "solid panipsilos"

c-------------
 98   Continue
c-------------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to see the whole surface"
      write (6,*) " 2 to see a particular elements"
      write (6,*) " 3 to compute the surface area and volume"
      write (6,*) " 4 to compute the surface area, volume"
      write (6,*) "   normal vector, and mean curvature"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c--------------------------
c prepare to print elements
c--------------------------

      if(menu.eq.1.or.menu.eq.2) then

       write (6,*)
       write (6,*) " Enter:"
       write (6,*)
       write (6,*) " 1 to print   six-node triangles"
       write (6,*) " 2 to print three-node triangles"
       write (6,*) " 0 to quit"
       write (6,*) " ---------"
       read  (5,*) Index

       if(Index.eq.0) Go to 98

      end If

c-----------------------
      if(menu.eq.1) then
c-----------------------

      if(Index.eq.1) then
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      else if(Index.eq.2) then
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end if

      Do k=1,nelm
       call printel (k,Index,color)      ! print in file "trgl6_hsph.net"
      End Do

      Do k=1,nelm
       call printel_stl (k,Index)  ! print in file "trgl6_hsph.stl"
      End Do

      Go to 98

c----------------------------
      else if(menu.eq.2) then
c----------------------------

      write (6,*)
      write (6,*) " Enter triangle number to see it "
      write (6,*) " and its neighbors"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) k

      if(k.eq.0) Go to 99

c---
c print the element and its neighbors
c---

      call printel (k,Index,color)
      call printel (nbe(k,1),Index,color)
      call printel (nbe(k,2),Index,color)
      call printel (nbe(k,3),Index,color)

      Go to 98

c-----------
      End If
c-----------

c------------------------------------
      if(menu.eq.3.or.menu.eq.4) then
c------------------------------------

c---
c Define surface integration
c triangle quadrature
c---

       write (6,*)
       write (6,*) " Will use the m-point gauss-triangle rule"
       write (6,*) " for integration over each triangle"
       write (6,*)
       write (6,*) " Please enter m"
       write (6,*)
       write (6,*) " choose from 1,3,4,6,7,9,12,13"
       write (6,*)
       write (6,*) " Enter 0 to quit"
       write (6,*) "----------------"
       read  (5,*) mint

       if(mint.eq.0) Go to 99

       call gauss_trgl
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

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

c-----------
      end if
c-----------


c-----------------------
      if(menu.eq.3) then
c-----------------------

c---
c Compute the surface area
c---

c---
c initialize
c---

      area = 0.0D0
      vlm  = 0.0D0

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

        Do i = 1,mint

          xi  = xiq(i)
          eta = etq(i)

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

          area = area +        hs*wq(i)
          vlm  = vlm  +  y*vny*hs*wq(i)

        End Do

      End Do

      area = 0.5D0*area      ! from the quadrature
      vlm  = 0.5D0*vlm       ! from the quadrature

      area_n = area/(pi2*req**2)        ! reduce
      vlm_n  = vlm /(pi2*req**3/3.0D0)    ! reduce

      write (6,*)
      write (6,110) area_n
      write (6,111) vlm_n
      write (6,*)

c----------------------------
      Else If(menu.eq.4) then
c----------------------------

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
     +  (nelm,npts,mint
     +  ,xmom,ymom,zmom
     +  ,area,vlm
     +  ,cx,cy,cz
     +  )

c---
c normalize and display
c---

      area_n = area/(pi2*req**2)        ! normalize
      vlm_n  = vlm /(pi2*req**3/3.0)    ! normalize

      write (6,*)
      write (6,110) area_n
      write (6,111) vlm_n
      write (6,*)

      write (6,*)
      write (6,*) " Element surface area and mean curvature"
      write (6,*)

      Do i=1,nelm
       write (6,100) i,arel(i),crvmel(i)
      End Do

c---------------------------
c Compute the mean curvature
c at the nodes
c
c Note that subroutine elm_geom 
c should be called before crvm6
c---------------------------

      call crvm6 (nelm,npts)

      write (6,*)
      write (6,*) " Node position, normal vector, mean curvature"
      write (6,*)

      Do i=1,npts
       write (6,102) i,p(i,1),p(i,2),p(i,3)
     +                ,vna(i,1),vna(i,2),vna(i,3)
     +                ,crvm(i)
      End Do

c-----------
      End If
c-----------


      Go to 98    ! return to the menu

c-------
c Ending
c-------

  99  Continue

      write (1,100) null
      close (1)

      write (8,*) "endsolid panipsilos"
      close (8)

c-----
c Done
c-----

  100 Format (1x,i4,10(1x,f12.5))
  101 Format (10(1x,f12.5))
  102 Format (1x,i4,10(1x,f8.5))

  110 Format (" Reduced surface Area :",F15.10)
  111 Format (" Reduced volume       :",F15.10)

      Stop
      End
