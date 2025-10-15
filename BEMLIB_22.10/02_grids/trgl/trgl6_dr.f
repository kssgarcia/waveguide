      program trgl6_dr

c=========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c---------------------------------------------------
c  Driver for discretizing a closed surface
c  into a network of 6-node curved triangles
c  based on the successive subdivision of
c  a regular octahedron (8 faces)
c  or icosahedron (20 faces)
c
c  This program generates an ordered grid
c  of six-node triangular elements
c  on the surface of an ellipsoid
c  and then performs various computations.
c
c  SYMBOLS:
c  -------
c
c  cxp, cyp, czp: coordinates of the center of the ellipsoid
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
c
c  crvt: curvature tensor at the nodes
c
c  rmat:   auxiliary matrix for holding the
c          modified Legendre function
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension     p(1026,3)
      Dimension    ne(1026,7)
      Dimension   vna(1026,3)
      Dimension  crvm(1026)
      Dimension  crvt(1026,3,3)
      Dimension color(1026)

      Dimension    u(1026,3)

      Dimension  dilt(1026)

      Dimension     fnc(1026),fnc1(1026),fnc2(1026)

      Dimension pxx(1026),pxy(1026),pxz(1026)
      Dimension           pyy(1026),pyz(1026)
      Dimension                     pzz(1026)

      Dimension       n(512,6),          nbe(512,3)
      Dimension   alpha(512),beta(512),gamma(512)
      Dimension    arel(512),xmom(512), ymom(512),zmom(512)
      Dimension  crvmel(512)

      Dimension xiq(20),etq(20),wq(20)

      Dimension rmat(0:32,0:32)   ! auxiliary matrix for the Leg fnc

      Dimension  v_aux(1026,3)    ! auxiliary vector
      Dimension    vtg(1026,3,3)  ! tangential gradient of a vector

      Dimension Gcm(1026,1026),Gmm(1026,1026)   ! conductivity and mass matrix

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/veloc0/u

      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna
      common/geo3/crvm
      common/geo6/crvmel
      common/geo9/xmom,ymom,zmom

      common/trq/xiq,etq,wq

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

      null    = 0
      Nfour   = 4
      Nseven  = 7

      oot = 1.0D0/3.0D0

c------------
c preferences
c------------

  76  Continue

      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to read input from file: trgl6.dat"
      write (6,*) " 2 to enter data"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ienrd

      if(Ienrd.eq.0) Go to 99

c-------------------------
      if(Ienrd.eq.1) then
c-------------------------

      open (9,file="trgl6.dat")

       read (9,*) Ishape
       read (9,*) Ndiv
       read (9,*) boa,coa
       read (9,*) req
       read (9,*) cxp,cyp,czp

      close (9)

c-----------------------------
      else if(Ienrd.eq.2) then
c-----------------------------

      write (6,*)
      write (6,*) " Will triangulate the surface of an ellipsoid"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for triangulation based on the octahedron"
      write (6,*) " 2 for triangulation based on the icosahedron"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ishape

      if(Ishape.eq.0) Go to 99

  91  write (6,*)
      write (6,*) " Will discretize the surface of an ellipsoid"
      write (6,*) " Please select the level of triangulation"
      write (6,*)
      write (6,*) " Enter 0, 1, 2, 3; 99 to quit"
      write (6,*) 
      write (6,*) " 0 will give an octahedron or icosahedron"
      write (6,*) " ----------------------------------------"
      read  (5,*) Ndiv 

      if(Ndiv.eq.99) Go to 99

      if(Ndiv.gt.3)  then
        write(6,*) 'Out of range: please try again.'
        Go to 91
      end if

      write (6,*)
      write (6,*) " The x, y, and z axes of the ellipsoid are"
      write (6,*) " denoted by a, b ,c"
      write (6,*)
      write (6,*) " Enter the axes ratios b/a and c/a"
      write (6,*) " ---------------------------------"
      read  (5,*) boa,coa

      write (6,*)
      write (6,*) " Enter the equivalent radius"
      write (6,*) " (defined with respect to the volume)" 
      write (6,*) " ------------------------------------"
      read  (5,*) req

      write (6,*)
      write (6,*) " Enter the coordinates of the"
      write (6,*) " center of the ellipsoid"
      write (6,*) " -----------------------"
      read  (5,*) cxp,cyp,czp

c---------
      Else
c---------

      write (6,*)
      write (6,*) " trgl: invalid choice"
      write (6,*)

      Go to 76

c-----------
      end if
c-----------

c----------------------------
c triangulate the unit sphere
c----------------------------

      if(Ishape.eq.1) call trgl6_octa (Ndiv,Npts,Nelm)
      if(Ishape.eq.2) call trgl6_icos (Ndiv,Npts,Nelm)
      
      write (6,*)
      write (6,*) " trgl: Number of nodes:    ",Npts
      write (6,*)
      write (6,*) " trgl: Number of elements: ",Nelm
      write (6,*)

c--------------------------------------
c expand to specified equivalent radius
c
c translate the center of the sphere
c to a specified position
c--------------------------------------

      scale = req/(boa*coa)**oot

      Do i=1,npts
        p(i,1) = scale*    p(i,1) + cxp
        p(i,2) = scale*boa*p(i,2) + cyp
        p(i,3) = scale*coa*p(i,3) + czp
      End Do

c--------------
c surface color
c--------------

      Do i=1,npts
       color(i) = p(i,1)
      End Do

c------------------
c open output files
c------------------

      open (1,file="trgl6.net")
      open (8,file="trgl6.stl")

      write (8,*) "solid deluge"

c-----------------------------

 98   Continue

      write (6,*) 
      write (6,*) "      MENU"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 to see the whole grid"
      write (6,*) " 2 to see a particular element"
      write (6,*) " 3 to compute the surface area and volume"
      write (6,*) " 4 to compute the surface area, volume,"
      write (6,*) "   normal vector, and mean curvature"
      write (6,*) 
      write (6,*) " 5 to compute the surface integral of a function"
      write (6,*) "   using a triangle quadrature"
      write (6,*) "51 to compute the surface integral of a function"
      write (6,*) "   using an integration rule"
      write (6,*) 
      write (6,*) " 6 to verify the orthonormality of"
      write (6,*) "   the modified Legendre functions"
      write (6,*) "   over the surface of a sphere"
      write (6,*) " 7 to compute the curvature tensor"
      write (6,*) " 8 to compute the surface divergence"
      write (6,*) "   of a velocity field (rate of strethcing)"
      write (6,*) " 9 for the conductivity and mass matrix"
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

      end if

c-----------------------
      if(menu.eq.1) then   ! print all triangles in file: trgl6.net
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
       call printel (k,Index,color) ! print in file "trgl6.net"
      End Do

      Do k=1,nelm
       call printel_stl (k,Index)  ! print in file "trgl6.stl"
      End Do

      Go to 98

c----------------------------
      else if(menu.eq.2) then
c----------------------------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " the label of a triangle to be seen"
      write (6,*) "     together with its neighbors"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) k

      if(k.eq.0) Go to 99

c---
c print the kth element
c and its 3 neighbors
c---

      call printel (k,Index,color)            ! will print in file "trgl6.net"
      call printel (nbe(k,1),Index,color)
      call printel (nbe(k,2),Index,color)
      call printel (nbe(k,3),Index,color)

      Go to 98

c-----------
      end if
c-----------

c---------------------------------------
c Perform miscellaneous computations
c based on the parametric representation
c---------------------------------------

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
     +  (p(i1,1),p(i1,2),p(i1,3)
     +  ,p(i2,1),p(i2,2),p(i2,3)
     +  ,p(i3,1),p(i3,2),p(i3,3)
     +  ,p(i4,1),p(i4,2),p(i4,3)
     +  ,p(i5,1),p(i5,2),p(i5,3)
     +  ,p(i6,1),p(i6,2),p(i6,3)
     +  ,alpha(k),beta(k),gamma(k)
     +  )

       End Do

c---------------------------
c select triangle quadrature
c for surface integration
c---------------------------

       if(menu.eq.3.or.menu.eq.4
     +             .or.menu.eq.5
     +             .or.menu.eq.6
     +             .or.menu.eq.7) then

       write (6,*)
       write (6,*) " Will use the m-point Gauss-triangle rule"
       write (6,*) "      for integration over each triangle"
       write (6,*)
       write (6,*) " Please enter m"
       write (6,*)
       write (6,*) " choose from 1,3,4,6,7,9,12,13"
       write (6,*)
       write (6,*) " 0 to quit"
       write (6,*) "----------"
       read  (5,*) mint

       if(mint.eq.0) Go to 99

       call gauss_trgl
     +
     +  (mint
     +  ,xiq,etq,wq
     +  )

       end if

c-----------------------
      if(menu.eq.3) then
c-----------------------

c---
c compute the surface area and volume
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
c loop over integration points
c---

        Do i=1,mint

          xi  = xiq(i)
          eta = etq(i)

c---
c interpolate the position vector
c             the normal vector
c             the surface metric
c---

          call interp_vn_hs
     +
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +
     +    ,alpha(k),beta(k),gamma(k)
     +    ,xi,eta
     +
     +    ,x,y,z
     +    ,vnx,vny,vnz
     +    ,hs
     +    )

          area = area +                     hs*wq(i)
          vlm  = vlm  + (x*vnx+y*vny+z*vnz)*hs*wq(i)/3.0D0

        End Do

      End Do

      area = 0.5D0 * area     ! factor 0.5 from the quadrature
      vlm  = 0.5D0 * vlm      ! factor 0.5 from the quadrature

      arean = area/(pi4*req**2)          ! normalize
      vlmn  = vlm /(pi4*req**3/3.0D0)    ! normalize

      write (6,*)
      write (6,110) arean
      write (6,111) vlmn
      write (6,*)

c----------------------------
      else if(menu.eq.4) then
c----------------------------

c----------------------------
c Compute:
c
c (a) the surface area of the individual elements.
c (b) the x, y, and z surface moments over each element.
c (c) the total surface area and volume.
c (d) the mean curvature of each element.
c (e) the average value of the normal vector at each node;
c     this is done by computing the normal vector at the 6 nodes
c     of each triangle, and then averaging the contributions.
c----------------------------

      call elm_geom
     +
     +   (nelm,npts
     +   ,mint
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c---
c normalize and display
c---

      arean = area/(pi4*req**2)        ! normalize
      vlmn  = vlm /(pi4*req**3/3.0D0)    ! normalize

      write (6,*)
      write (6,110) arean
      write (6,111) vlmn
      write (6,112) cx,cy,cz
      write (6,*)

      write (6,*) 
      write (6,*) " Element surface area and element mean curvature"
      write (6,*) 

      Do i=1,nelm
       write (6,100) i,arel(i),crvmel(i)
      End Do

c------
c Compute the mean curvature at the nodes
c
c Note that subroutine crvm6
c should be called after elm_geom
c------

      call crvm6 (nelm,npts)

      write (6,*) 
      write (6,*) " Node position, normal vector, mean curvature"
      write (6,*) 

      Do i=1,npts
       write (6,102) i,p(i,1),p(i,2),p(i,3)
     +                ,vna(i,1),vna(i,2),vna(i,3)
     +                ,crvm(i)
      End Do

c------------------------------------------
      else if(menu.eq.5.or.menu.eq.51) then
c------------------------------------------

c----------------------------
c Compute the surface integral
c of a function
c----------------------------

 93   Continue

      write (6,*) 
      write (6,*) " MENU OF FUNCTIONS"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for f = 1.0"
      write (6,*) " 2 for f = x**2"
      write (6,*) " 3 for f = x**3"
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read  (5,*) my_choice

      If(my_choice.eq.0) Go to 98

c---
c evaluate the integrand at the nodes
c---
     
      Do i=1,npts
        If(my_choice.eq.1) fnc(i) = 1.0D0
        If(my_choice.eq.2) fnc(i) = p(i,1)**2
        If(my_choice.eq.3) fnc(i) = p(i,1)**3
      End Do

c---
c perform the integration
c---

c---
      if(menu.eq.5) then    ! integration quadrature

       call srf_int_3d
     +
     +    (nelm,npts
     +    ,mint
     +    ,fnc
     +    ,olok
     +    )

      else if(menu.eq.51) then ! integration rule 

       write (6,*) " Enter the order of the integration rule"
       write (6,*) " 0 to quit"
       write (6,*) "----------"
       read  (5,*) mint1

       if(mint1.eq.0) Go to 99

       call srf_int_3d_ir
     +
     +    (nelm,npts
     +    ,mint1
     +    ,fnc
     +    ,olok
     +    )

      end if
c---

      write (6,*)
      write (6,456) olok

      Go to 93

c----------------------------
      else if(menu.eq.6) then
c----------------------------

c-------------------------------------------
c Verify the orthonormality of the modified
c associated Legendre functions over the
c surface of a sphere
c-------------------------------------------

      write (6,*) 
      write (6,*) " Please enter the maximum order j_max"
      write (6,*) "         must be 0,1,3,..."
      write (6,*) " Enter 98 to quit"
      write (6,*) 
      write (6,*) " Legendre functions P_j^m will be computed"
      write (6,*) " for j=0,...,j_max and m=0,...,j"
      write (6,*) 
      write (6,*) " Maximum j_max is 32"
      write (6,*) " -------------------------------"
      read  (5,*) j_max

      If(j_max.eq.98) Go to 98 ! return to the main menu

  84  write (6,*) 
      write (6,*) " Testing orthonormality"
      write (6,*) 
      write (6,*) " Please enter two pairs (j1, m1), (j2,m2)"
      write (6,*) " j1, j2 =0,...,",j_max
      write (6,*) 
      write (6,*) " Will stop if j1=98 or j2=98 "
      write (6,*) " ----------------------------"
      read  (5,*) j1,m1,j2,m2

      If(j1.eq.98) Go to 98 ! return to the main menu
      If(j2.eq.98) Go to 98 ! return to the main menu

c---
c trap
c---

      if((j1.lt.0).or.(j2.lt.0.).or.(j1.gt.j_max)
     +                          .or.(j2.gt.j_max)
     +  ) then

        write (6,*)
        write (6,*) "Not possible; please try again"
        write (6,*)
        Go to 84

      end if
     
c---
c evaluate the integrand at the nodes
c---

      Do i=1,npts

        rdist = sqrt( p(i,1)**2+p(i,2)**2+p(i,3)**2 )
        costh = p(i,1)/rdist

        call leg_fnc_ortho 
     +
     +     (costh
     +     ,j_max
     +     ,rmat
     +     )

        sigma  = sqrt( p(i,2)**2+p(i,3)**2 )

        if(sigma.lt.0.0000000001) then
         phi = 0.0D0
        else
         cosphi = p(i,2)/sigma
         if(cosphi.gt. 0.999999999) cosphi = 0.999999999
         if(cosphi.lt.-0.999999999) cosphi =-0.999999999
         phi = acos(cosphi)
         if(p(i,3).lt.0) phi = pi2-phi
        end if

        first_c = rmat(j1,m1)*cos(m1*phi)
        first_s = rmat(j1,m1)*sin(m1*phi)
        secon_c = rmat(j2,m2)*cos(m2*phi)
        secon_s = rmat(j2,m2)*sin(m2*phi)

        fnc1(i) = first_c*secon_c + first_s*secon_s
        fnc2(i) = first_c*secon_s - first_s*secon_c

      End Do

c---
c perform the surface integration
c---

      call srf_int_3d
     +
     +    (nelm,npts,mint
     +    ,fnc1
     +    ,olokr
     +    )

      call srf_int_3d
     +
     +    (nelm,npts,mint
     +    ,fnc2
     +    ,oloki
     +    )

      write (6,*)
      write (6,*) " Real and Imaginary parts " 
      write (6,*)
      write (6,457) olokr,oloki

      Go to 84

c============================
      else if(menu.eq.7) then
c============================

c-----------------------------
c compute the curvature tensor
c-----------------------------

      call elm_geom
     +
     +  (nelm,npts
     +  ,mint
     +  ,area,vlm
     +  ,cx,cy,cz
     +  )

c---
c compute the curvature tensor
c---

      Do i=1,npts
       v_aux(i,1) = vna(i,1)     ! auxiliary vector
       v_aux(i,2) = vna(i,2)
       v_aux(i,3) = vna(i,3)
      End Do

c------
c Note that elm_geom should be called
c called before vtg_3d
c------

      call vtg_3d
     +
     +  (nelm
     +  ,npts
     +  ,v_aux
     +  ,crvt
     +  )

c---------------------------
c compute 2km = Trace{crvt}
c (twice the mean curvature)
c---------------------------

      Do i=1,npts

       trace = crvt(i,1,1)+crvt(i,2,2)+crvt(i,3,3)

       crv_mean = 0.50D0*trace

c------------
c compute B.n
c------------

       prjx = 0.0D0
       prjy = 0.0D0
       prjz = 0.0D0

       Do k=1,3
        prjx = prjx + crvt(i,1,k) * vna(i,k)
        prjy = prjy + crvt(i,2,k) * vna(i,k)
        prjz = prjz + crvt(i,3,k) * vna(i,k)
       End Do

c------------
c compute n.B
c------------

       prjx1 = 0.0D0
       prjy1 = 0.0D0
       prjz1 = 0.0D0

       Do k=1,3
        prjx1 = prjx1 + crvt(i,k,1) * vna(i,k)
        prjy1 = prjy1 + crvt(i,k,2) * vna(i,k)
        prjz1 = prjz1 + crvt(i,k,3) * vna(i,k)
       End Do

c-------------
c compute I-nn
c-------------

       pxx(i) = 1.0D0 - vna(i,1)**2
       pxy(i) =       - vna(i,1)*vna(i,2)
       pxz(i) =       - vna(i,1)*vna(i,3)
       pyy(i) = 1.0D0 - vna(i,2)**2
       pyz(i) =       - vna(i,2)*vna(i,3)
       pzz(i) = 1.0D0 - vna(i,3)**2

c---------
c printing
c---------

       write (6,100) i,crv_mean

c      write (6,103) crvt(i,1,1),crvt(i,1,2)
c    +              ,crvt(i,1,3),prjx,prjx1
c      write (6,103) crvt(i,2,1),crvt(i,2,2)
c    +              ,crvt(i,2,3),prjy,prjy1
c      write (6,103) crvt(i,3,1),crvt(i,3,2)
c    +              ,crvt(i,3,3),prjz,prjz1

c      write (6,103) pxx(i),pxy(i),pxz(i)
c      write (6,103) pxy(i),pyy(i),pyz(i)
c      write (6,103) pxz(i),pyz(i),pzz(i)

      End Do

c============================
      else if(menu.eq.8) then
c============================

c---
c specify the velocity
c---

      Do i=1,Npts
c      u(i,1) = p(i,2)
c      u(i,2) = 0.0D0
c      u(i,3) = 0.0D0
       u(i,1) = p(i,1)
       u(i,2) = p(i,2)
       u(i,3) = p(i,3)
      End Do

      call dilt_3d
     +
     +  (Nelm
     +  ,Npts
     +  ,dilt
     +  )

      write (6,*) 
      write (6,*) " trgl6: rate of surface dilatation:"
      write (6,*) 

      Do i=1,Npts
       write (6,100), i,p(i,1),dilt(i)
      End Do

c============================
      else if(menu.eq.9) then
c============================

      call cmm
     +
     +   (Npts,Nelm
     +   ,mint
     +   ,Gcm,Gmm
     +   )

      Do i=1,Npts
        Do j=1,Npts
         write (6,104) i,j,Gcm(i,j),Gmm(i,j)
        End Do
      End Do

c============
      End If
c============

      Go to 98    ! return to the menu

c-------
c ending
c-------

  99  Continue

      write (1,100) null
      close (1)

      write (8,*) "endsolid deluge"
      close (8)

c-----
c Done
c-----

  100 Format (1x,i4,10(1x,f12.5))
  101 Format (10(1x,f12.5))
  102 Format (1x,i4,10(1x,f8.5))
  103 Format (10(1x,f8.5))
  104 Format (1x,i4,1x,i4,10(1x,f10.5))

  110 Format (" Surface Area :",F15.10)
  111 Format (" Volume       :",F15.10)
  112 Format (" Centroid: ",F15.10,1x,F15.10,1x,F15.10)

  456 Format ("Integral = ",f15.10)
  457 Format ("Integrals = ",f15.10,1x,f15.10)

      Stop
      End
