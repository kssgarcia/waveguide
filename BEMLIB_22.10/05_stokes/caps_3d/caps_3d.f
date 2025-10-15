      program caps_3d

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-----------------------------------------------
c Dynamic simulation of the deformation of a liquid
c capsule subject to an imposed simple shear flow
c
c Three types of flow:
c
c (a) Flow in an infinite domain
c (b) Flow a semi-infinite domain bounded by a plane wall
c (c) Triply periodic flow
c
c In this implementation, the surface tension is constant.
c
c Surfactant concentration is assigned to the nodes,
c but is dynamically inactive; that is, the surfactant
c concentration remains constant in time.
c
c The interface develops
c elastic tensions and bending moments
c
c SYMBOLS:
c --------
c
c  npts	   total number of points
c  nelm	   total number of elements
c
c  p(i,j)   coordinates of nodes i (j=1,2,3)
c
c  ne(k,j)  ne(k,1) is the number of elements adjacent to point k
c           ne(k,2), ... are the elements numbers, j = 2, ..., 7
c           For this triangulation, ne(k,j) is up to six
c
c  n(k,i)    connectivity table: points for element k, i = 1,...,6
c
c  nbe(k,j)  the three neighboring elements of element k (j=1,2,3)
c
c  vna       average value of the normal vector at nodes
c  u         velocity of points
c  c(i)      surfactant concentration at node i   
c  srtn(i)   isotropic surface tension at node i   
c
c  ps                          save p
c  Usave,Vsave,Wsave           save U,V,W
c
c    arel(k)    surface area of element k
c  crvmel(k)    average value of the mean curvature
c               of element k
c    xmom(k)    x-moment of element k
c    ymom(k)    y-moment of element k
c    zmom(k)    z-moment of element k
c
c  jxy          number of sequential nodes in the xy plane
c
c  dxy          Taylor deformation parameter in the xy plane
c  thmax        inclination of maximum axis
c  thmin        inclination of minimum axis
c
c  vol         capsule volume
c  ars         capsule surface area 
c  crx         x position of centroid   
c  cry         y position of centroid   
c  crz         z position of centroid   
c
c  zz, ww         base points and weights of Gauss-Legenrde integr
c  xiq, etq, wq   base points and weights of Gauss-triangle integr
c
c  wall    wall located at y = wall
c
c Elst:      modules of elasticity
c Elstb:     bending modules of elasticity
c crvmr:     mean curvature of the reference shape
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension     p(1026,3)
      Dimension    ne(1026,7)
      Dimension   vna(1026,3)
      Dimension     u(1026,3)
      Dimension     c(1026),srtn(1026)
      Dimension Umove(1026),Vmove(1026),Wmove(1026)

      Dimension   pr(1026,3)   ! reference position
      Dimension vnar(1026,3)   ! reference normal vector

      Dimension ps(1026,3),Usave(1026),Vsave(1026),Wsave(1026)

      Dimension nvel(1026),lxy(1026,2)

      Dimension crvm(1026)

      Dimension      n(512,6), nbe(512,3)
      Dimension  alpha(512),  beta(512), gamma(512)
      Dimension   arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension crvmel(512)

      Dimension jxy(100)
      Dimension x_xy(100),y_xy(100),ux_xy(100),uy_xy(100)

      Dimension time(2000),dxy(2000),thmax(2000),thmin(1000)
      Dimension  ars(2000),vol(2000)
      Dimension  crx(2000),cry(2000),  crz(2000)
      Dimension  sxy(2000),sd1(2000),  sd2(2000)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna
      common/geo3/crvm
      common/geo6/crvmel
      common/geo8/nxy,jxy
      common/geo9/xmom,ymom,zmom

      common/viscr/vs1,vs2,vsrt,vsrtp,vsrtm,vsf,vsk
      common/visci/Ivs

      common/var/shrt,wall

      common/veloc1/nvelt,nvel
      common/veloc2/nvelr,lxy
      common/tension/srtn

      common/elten1/pr,vnar

c---
c for doubly-periodic flow
c---

      common/aaaa/a11,a12,a13,a21,a22,a23,a31,a32,a33
      common/bbbb/b11,b12,b13,b21,b22,b23,b31,b32,b33
      common/ewew/ew,tau
      common/mmmm/Max1,Max2

c---
c various
c---

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c---
c constants
c---

      pi = 3.14159 265358D0

      piq = 0.25D0*pi
      pih = 0.50D0*pi
      pi2 = 2.00D0*pi
      pi4 = 4.00D0*pi
      pi6 = 6.00D0*pi
      pi8 = 8.00D0*pi

      Null = 0
      None = 1
      Nfour = 4
      Nseven = 7

      oot = 1.0D0/3.0D0

c------
c input
c------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to read them from file: caps_3d.dat'
      write (6,*) ' 2 to type the data into the keyboard'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*)   Ienrd

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (2,file="caps_3d.dat")

      read (2,*) Ioctaicos
      read (2,*) ndiv 
      read (2,*) 
      read (2,*) Ishape
      read (2,*) req
      read (2,*) boa,coa
      read (2,*) cx,cy,cz
      read (2,*) phi1,phi2,phi3
      read (2,*)
      read (2,*) mint
      read (2,*) NGL
      read (2,*)
      read (2,*) Iflow
      read (2,*)
      read (2,*) wall
      read (2,*)
      read (2,*) a11,a12,a13
      read (2,*) a21,a22,a23
      read (2,*) a31,a32,a33
      read (2,*) Max1,Max2
      read (2,*)
      read (2,*) shrt
      read (2,*)
      read (2,*) vs1
      read (2,*) vs2
      read (2,*)
      read (2,*) srft
      read (2,*) 
      read (2,*) Elst
      read (2,*) Elstb
      read (2,*) crvmr
      read (2,*)
      read (2,*) nter
      read (2,*) tol
      read (2,*) Idfl
      read (2,*)
      read (2,*) Iread
      read (2,*)
      read (2,*) norm
      read (2,*) Isym_xy
      read (2,*)
      read (2,*) IRK
      read (2,*) Dt
      read (2,*)
      read (2,*) Nprint_xy
      read (2,*) Nprint_xyz
      read (2,*)
      read (2,*) move
      read (2,*)

c------------------------
      else if(Ienrd.eq.2) then
c------------------------

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to discretize an octahedron"
      write (6,*) " 2 to discretize an icosahedron"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ioctaicos

      if(Ioctaicos.eq.0) stop

 88   write (6,*)
      write (6,*)
      write (6,*) " Enter the level of triangulation: 0, 1, 2, 3"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) ndiv 

      if(ndiv.eq.99) Go to 99

      if(ndiv.gt.3) then
        write (6,*) 'Level is too high, should be less than 4'
        Go to 88
      end if

      write (6,*)
      write (6,*) " Choose the unstressed shape:"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for an ellipsoid"
      write (6,*) " 2 for an RBC"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ishape

      if(Ishape.eq.99) Go to 99

c-----
      if(Ishape.eq.1) then  !  SPHEROID
c-----

       write (6,*)
       write (6,*) " Enter the equivalent radius"
       write (6,*) " ---------------------------"
       read  (5,*) req

       write (6,*)
       write (6,*) " The capsule has an ellipsoidal initial shape"
       write (6,*) " with x,y,z semi-axis: a,b,c"
       write (6,*)
       write (6,*) " Enter the axes ratios b/a and c/a"
       write (6,*) " ---------------------------------"
       read  (5,*) boa,coa

c-----
      end if
c-----

      write (6,*)
      write (6,*) " Enter coordinates of the capsule center"
      write (6,*) " ---------------------------------------"
      read  (5,*) cxp,cpy,czp

      write (6,*)
      write (6,*) " Enter the three rotation angles about the"
      write (6,*) "       x,y,z axes, in multiples of pi"
      write (6,*) " -----------------------------------------"
      read  (5,*) phi1,phi2,phi3

      write (6,*)
      write (6,*) " Non-singular integration over each triangle"
      write (6,*)
      write (6,*) " Will use the m-point rule."
      write (6,*) 
      write (6,*) " Please enter m"
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) mint

      if(mint.eq.0) Go to 99

      write (6,*)
      write (6,*) " Singular integration over each triangle"
      write (6,*)
      write (6,*) " Will use the m-point Gauss-Legendre quadrature"
      write (6,*) 
      write (6,*) " Please enter m"
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*) " 0 to quit"
      write (6,*) " ----------"
      read  (5,*) NGL

      if(NGL.eq.0) Go to 99

      write (6,*)
      write (6,*) " Choose the flow"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for flow in infinite space"
      write (6,*) " 2 for flow bounded by a plane wall"
      write (6,*) " 3 for triply-periodic flow"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Iflow

      if(Iflow.eq.0) Go to 99

c----
c wall bounded flow
c----

      if(Iflow.eq.2) then

       write (6,*)
       write (6,*) " Wall is located at y = wall"
       write (6,*)
       write (6,*) " Please enter: wall"
       write (6,*) " ------------------"
       read  (5,*) wall

      end if

c----
c triply-periodic flow
c----

      if(Iflow.eq.3) then

       write (6,*) " Triply-periodic flow"
       write (6,*)
       write (6,*) ' Enter the x, y, and z coordinates'
       write (6,*) '       of the first lattice vector'
       write (6,*) ' ---------------------------------'
       read  (5,*)  a11,a12,a13

       write (6,*) ' Enter the x, y, and z coordinates'
       write (6,*) '       of the second lattice vector'
       write (6,*) ' ---------------------------------'
       read  (5,*)  a21,a22,a23

       write (6,*) ' Enter the x, y, and z coordinates'
       write (6,*) '       of the third lattice vector'
       write (6,*) ' ---------------------------------'
       read  (5,*)  a31,a32,a33

       write (6,*)
       write (6,*) ' Enter the truncation limits Max1 and MAx2'
       write (6,*) '       for the Ewald summation of the'
       write (6,*) '       doubly-periodic Green function'
       write (6,*) '       in real and reciprocal space'
       write (6,*) ' ---------------------------------------'
       read  (5,*)   Max1,Max2

      end if

      write (6,*)
      write (6,*) " Simple shear flow"
      write (6,*)
      write (6,*) " Enter the shear rate" 
      write (6,*) " --------------------"
      read  (5,*) shrt

      write (6,*)
      write (6,*) " Enter the viscosity of the ambient fluid"
      write (6,*) " ----------------------------------------"
      read  (5,*) vs1

      write (6,*)
      write (6,*) " Enter the viscosity of the drop"
      write (6,*) " -------------------------------"
      read  (5,*) vs2

      write (6,*)
      write (6,*) " Enter the surface tension"
      write (6,*) " -------------------------"
      read  (5,*) srft

      write (6,*)
      write (6,*) " Enter the modulus of elastiticity"
      write (6,*) " ---------------------------------"
      read  (5,*) Elst

      write (6,*)
      write (6,*) " Enter the bending modulus of elastiticity"
      write (6,*) " -----------------------------------------"
      read  (5,*) Elstb

      write (6,*)
      write (6,*) " Enter the mean curvature of the reference state"
      write (6,*) " -----------------------------------------------"
      read  (5,*) crvmr

      write (6,*)
      write (6,*)
      write (6,*) " Enter the maximum number of iterations "
      write (6,*) "       in solving the integral equation"
      write (6,*) " --------------------------------------"
      write (6,*)
      write (6,*)
      write (6,*) " Enter the maximum number of iterations "
      write (6,*) "       in solving the integral equation"
      write (6,*) " --------------------------------------"
      read  (5,*) Nter

      write (6,*)
      write (6,*) " Enter the error tolerance "
      write (6,*) "       in solving the integral equation"
      write (6,*) " --------------------------------------"
      read  (5,*) tol

      write (6,*)
      write (6,*) " Enable deflation ? "
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for deflation of one eigenvalue "
      write (6,*) " 0 for no"
      write (6,*) " --------"
      read  (5,*) Idfl

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 to generate the initial shape"
      write (6,*) " 1 to read from restart file: caps_3d.inp "
      write (6,*) " -----------------------------------------"
      read  (5,*) Iread

      write (6,*)
      write (6,*) " Normalize volume after each step ?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      read  (5,*) norm

      write (6,*)
      write (6,*) " Exploit symmetry wr to the xy plane ?"
      write (6,*)
      write (6,*) " Enter 0 for no, 1 for yes"
      write (6,*) " -------------------------"
      read  (5,*) Isym_xy

      write (6,*)
      write (6,*) " Time integration"
      write (6,*)
      write (6,*) " Enter 1 for the Euler explicit method"
      write (6,*) "       2 for the RK2 method "
      write (6,*) " ---------------------------"
      read  (5,*) IRK

      write (6,*)
      write (6,*) " Enter the time step Dt"
      write (6,*) " ----------------------"
      read  (5,*) Dt

      write (6,*)
      write (6,*) " Will record a profile in the xy plane"
      write (6,*) " after N steps; please enter N"
      write (6,*) " -----------------------------"
      read  (5,*) Nprint_xy

      write (6,*)
      write (6,*) " Will record a 3D profile"
      write (6,*) " after N steps; please enter N"
      write (6,*) " -----------------------------"
      read  (5,*) Nprint_xyz

      write (6,*)
      write (6,*) ' Select node advancement method'
      write (6,*)
      write (6,*) ' Enter 0 to move with the total velocity '
      write (6,*) '       1 to move with the normal velocity '
      write (6,*) ' -----------------------------------------'
      read  (5,*) move

c-----------
      end if               ! End of reading parameters
c-----------

c     write (6,*) ndiv 
c     write (6,*) Ishape
c     write (6,*) req
c     write (6,*) boa,coa
c     write (6,*) cx,cy,cz
c     write (6,*) phi1,phi2,phi3
c     write (6,*) mint
c     write (6,*) NGL
c     write (6,*) Iflow
c     write (6,*) wall
c     write (6,*) shrt
c     write (6,*) vs1
c     write (6,*) vs2
c     write (6,*) srft
c     write (6,*) nter
c     write (6,*) tol
c     write (6,*) Idfl
c     write (6,*) Iread
c     write (6,*) norm
c     write (6,*) Isym_xy
c     write (6,*) IRK
c     write (6,*) Dt
c     write (6,*) Nstep
c     write (6,*) Nprint_xy
c     write (6,*) Nprint_xyz
c     write (6,*) move

c---------------------------------------
c preparations, adjustments, definitions
c---------------------------------------

c---------------------------
ctTriply-periodic shear flow
c---------------------------

      if(Iflow.eq.3) then

        a11h  = 0.5*a11    ! used for reseting
        a11hm = - a11h     ! the second lattice vector

      end if

c---
c set viscosity ratio index, etc
c---

      vsrt = vs2/vs1   ! viscosity ratio

c---
c viscosity ratio not equal to unity
c---

      Ivs = 0

      if(abs(vsrt-1.0).gt.0.0000001) then
       Ivs = 1
       vsrtm = 1.0D0-vsrt
       vsrtp = 1.0D0+vsrt
       vsf   = 2.0D0/vsrtp
       vsk   = vsrtm/vsrtp
      end if

c---
c infinite simple shear flow:
c place capsule center at the origin
c---

      if(Iflow.eq.1.or.Iflow.eq.3) then
        cxp = 0.0
        cyp = 0.0
        czp = 0.0
      end if

c-------------------------------
c Define integration quadratures
c-------------------------------

      call gauss_leg (NGL,zz,ww)

      call gauss_trgl 
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

c----------------------------
c triangulate the unit sphere
c
c Run even at restart
c to generate the connectivity table
c-----------------------------------

      if(Ioctaicos.eq.1) then

       call trgl6_octa
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      else

      call trgl6_icos
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      end if

      write (6,*)
      write (6,*) "Number of points:   ",npts
      write (6,*)
      write (6,*) "Number of elements: ",nelm
      write (6,*)

      nelmh = nelm/2
      nelm2 = 2*nelm

c----------------
c reference state 
c----------------

      time(1) = 0.0D0

c---------------------------
      if(Ishape.eq.1) then    ! SPHEROID
c---------------------------

c-------
c expand to specified shape
c and equivalent radius
c
c rotate by the angles: phi1,phi2,phi3
c translate center to specified position
c-------

      scale = req/(boa*coa)**oot

      Do i=1,npts
        p(i,1) = scale*    p(i,1)
        p(i,2) = scale*boa*p(i,2)
        p(i,3) = scale*coa*p(i,3)
      End do

      phi1 = phi1*pi
      phi2 = phi2*pi
      phi3 = phi3*pi

      cs = cos(phi1)
      sn = sin(phi1)

      Do i=1,npts
       tmpx = p(i,1)
       tmpy = cs*p(i,2)+sn*p(i,3)
       tmpz =-sn*p(i,2)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = cos(phi2)
      sn = sin(phi2)

      Do i=1,npts
       tmpx = cs*p(i,1)-sn*p(i,3)
       tmpy = p(i,2)
       tmpz = sn*p(i,1)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = cos(phi3)
      sn = sin(phi3)

      Do i=1,npts
       tmpx = cs*p(i,1)-sn*p(i,2)
       tmpy = sn*p(i,1)+cs*p(i,2)
       tmpz = p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      Do i=1,npts
        p(i,1) = p(i,1) + cxp
        p(i,2) = p(i,2) + cyp
        p(i,3) = p(i,3) + czp
      End Do

c---------
      else       ! RBC
c---------

c-------
c Deform unit sphere to a biconcave disk
c The disk is axially symmetric about the y axis
c--------
 
      Do i=1,npts
       ss = p(i,1)**2 + p(i,3)**2
       sq = ss**2
       p(i,2) = 0.5*p(i,2)*(0.207+2.003*ss-1.123*sq)
c      write (6,100) i,p(i,1),p(i,2),p(i,3),ss
      End Do

c-----
c expand to equivalent radius
c-----

      scale = 1.3858189*req

      Do i=1,npts
        p(i,1) = scale*p(i,1)
        p(i,2) = scale*p(i,2)
        p(i,3) = scale*p(i,3)
      End Do

c-------
c rotate by angles: phi1,phi2,phi3
c-------

      phi1 = phi1*pi
      phi2 = phi2*pi
      phi3 = phi3*pi

      cs = cos(phi1)
      sn = sin(phi1)

      Do i=1,npts
       tmpx = p(i,1)
       tmpy = cs*p(i,2)+sn*p(i,3)
       tmpz =-sn*p(i,2)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = cos(phi2)
      sn = sin(phi2)

      Do i=1,npts
       tmpx = cs*p(i,1)-sn*p(i,3)
       tmpy = p(i,2)
       tmpz = sn*p(i,1)+cs*p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      cs = cos(phi3)
      sn = sin(phi3)

      Do i=1,npts
       tmpx = cs*p(i,1)-sn*p(i,2)
       tmpy = sn*p(i,1)+cs*p(i,2)
       tmpz = p(i,3)
       p(i,1) = tmpx
       p(i,2) = tmpy
       p(i,3) = tmpz
      End Do

      Do i=1,npts
        p(i,1) = p(i,1) + cxp
        p(i,2) = p(i,2) + cyp
        p(i,3) = p(i,3) + czp
      End Do

c-----------
      end if      ! end of shape options
c-----------

c---
c unscale to record
c---

      phi1 = phi1/pi
      phi2 = phi2/pi
      phi3 = phi3/pi

c--------------------------------
c set the reference node position
c--------------------------------

      Do i=1,npts
       pr(i,1) = p(i,1)
       pr(i,2) = p(i,2)
       pr(i,3) = p(i,3)
      End Do

c------------------------------------
c compute the reference normal vector
c------------------------------------

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

        call elm_geom 
     +
     +    (nelm,npts,mint
     +    ,area,vlm
     +    ,cx,cy,cz
     +    )

      Do i=1,npts
       vnar(i,1) = vna(i,1)
       vnar(i,2) = vna(i,2)
       vnar(i,3) = vna(i,3)
      End Do

c-----------------------------------------
c Assign surfactant concentration to nodes
c
c Compute surface tension
c
c In this implementation, the surfactnat
c concentration remains constant in time
c-----------------------------------------

      Do i=1,npts
       c(i) = 1.0D0
       srtn(i) = srft
      End Do

c------------------------
c Initialize the velocity
c------------------------

      Do i=1,npts
       u(i,1) = 0.0D0
       u(i,2) = 0.0D0
       u(i,3) = 0.0D0
      End Do

c--------------------
      if(Iread.eq.1) then  ! read initial shape from file: caps_3d.inp
c--------------------

c-------------------------------------
c read coordinates of nodes and
c surfactant concentration from file: caps_3d.inp
c
c assign surface tension at nodes
c-------------------------------------

       open (9,file="caps_3d.inp")

       if(Iflow.eq.1.or.Iflow.eq.2) then
         read (9,*) npts,time(1)
       else if(Iflow.eq.3) then
         read (9,*) npts,time(1),a21
       end if

       Do i=1,npts
         read (9,*) idle,p(i,1),p(i,2),p(i,3)
         srtn(i) = srft
       End Do

       close (9)

c-----------
      end if
c-----------

c----
c compute the surface centroid
c---

      if(Iflow.eq.2) then

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

        call elm_geom 
     +
     +    (nelm,npts
     +    ,mint
     +    ,area,vlm
     +    ,cx,cy,cz
     +    )

      end if

c-----------------------------------------------
c identify points in the z=cz plane and order them
c
c nxy is number of points in the xy plane
c the vector jxy contains consecutive points
c-----------------------------------------------

      if(Isym_xy.eq.1) then

       call xy_slice 
     +
     +   (npts
     +   ,cx,cy,cz
     +   ,nxy
     +   ,jxy
     +   )

       Do i=1,nxy
        write (6,*) i,jxy(i)
       End Do

      end if

c------------------------------------------
c Prepare to exploit symmetry with respect
c to the xy plane
c
c nvelt:  number of nodes where the velocity 
c                           will be computed
c nvelr:  number of nodes where the velocity 
c                          will be reflected 
c
c lxy(i,1): number of reflected node i
c lxy(i,2): its image
c------------------------------------------

c---
      if(Isym_xy.eq.0) then   ! will compute velocity at all nodes
c---

       nvelt = npts

       Do i=1,npts
        nvel(i) = i
       End Do

c---
      else       ! will compute the velocity at roughly half nodes
c---

      icount1 = 0
      icount2 = 0

      Do 66 i=1,npts

        if(p(i,3).gt.-0.0001) then
          icount1 = icount1+1
          nvel(icount1) = i
        else
          icount2 = icount2+1
          lxy(icount2,1) = i
          Do j=1,npts
           test1 = abs(p(i,1)-p(j,1))
           test2 = abs(p(i,2)-p(j,2))
           test3 = abs(p(i,3)+p(j,3))
           if(    test1.le.0.0001
     +       .and.test2.le.0.0001
     +       .and.test3.le.0.0001
     +       ) then
             lxy(icount2,2) = j
             Go to 66
           end if
          End Do
        End If

  66  Continue

      nvelt = icount1
      nvelr = icount2

c---
      end if
c---

c------------
c born to run
c------------

      open (3,file="caps_3d.xy")
      open (4,file="caps_3d.diag")
      open (8,file="caps_3d.xyz")

      Kstep = 1           ! time step counter

      Iprint_xy  = Nprint_xy       ! printing counter
      Iprint_xyz = Nprint_xyz      ! printing counter

c-------------------------------------
c generate a Matlab visualization file
c-------------------------------------

      open (1,file="caps_3d.net")

      Index = 1  ! 6-node triangles
      Index = 2  ! 3-node triangles

      if(Index.eq.1) then          ! 6-node triangles
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      else if(Index.eq.2) then     ! 3-node triangles
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end if

      write (1,*) Iflow
      write (1,*) wall

c-.-.-.-.-.-.-.-.-.-.-.-.-.-.
c  time stepping begins here
c-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 90   Continue

      if(Ienrd.eq.0) then
        write (6,*)
        write (6,*) 'Enter number of steps before pausing'
        write (6,*) '0 to quit'
        write (6,*) '---------'
        read  (5,*)  Nstep
      else
        read  (2,*)  Nstep
      end if

      if(Nstep.eq.0) Go to 99

c-----------
c initialize
c-----------

      Istep  = 1      ! batch step counter

  97  Continue

      write (6,*)
      write (6,*) "--------------------"
      write (6,*)
      write (6,105) Istep,Nstep,time(Kstep)

      Ipass = 0      ! Ipass will be used for volume normalization

  91  Continue

      Ipass = Ipass + 1 

c----------------------------------------
c compute coefficients alpha, beta, gamma
c for quadratic xi-eta mapping
c over each element
c---------------------------------------

      Do k=1,nelm

        i1 = n(k,1)
        i2 = n(k,2)
        i3 = n(k,3)
        i4 = n(k,4)
        i5 = n(k,5)
        i6 = n(k,6)

        call abc 
     +
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +    ,alpha(k),beta(k),gamma(k)
     +    )

      End Do

c------------------------------------------------
c compute:
c
c     surface area of the individual elements
c     x, y, and z surface moments over each element
c     total surface area and volume
c     mean curvature over each element
c     average normal vector at the nodes
c------------------------------------------------

      call elm_geom 
     +
     +  (nelm,npts,mint
     +  ,area,vlm
     +  ,cx,cy,cz
     +  )

c---------------
c     write (6,*) 
c     write (6,*) " Curvature of elements"
c     write (6,*) 
c     Do k = 1,nelm
c       write (6,100) k,crvmel(k)
c     End Do
c---------------

c--------------------------
c normalize capsule volume?
c--------------------------

      if(Norm.eq.1) then

       if(Ipass.eq.1) then 

         cf = (pi4*req**3/(3.0*vlm))**oot

         Do i=1,npts
           p(i,1) = p(i,1)*cf
           p(i,2) = p(i,2)*cf
           p(i,3) = p(i,3)*cf
         End Do

         Go to 91

       end if

      end if

c------------------
c recording session
c------------------

      ars(Kstep) = area/(pi4*req**2)
      vol(Kstep) = vlm/ (pi4*req**3/3.0)

      crx(Kstep) = cx
      cry(Kstep) = cy
      crz(Kstep) = cz

      write (6,*)
      write (6,110) ars(istep)
      write (6,111) vol(istep)
      write (6,112) cx,cy,cz
      write (6,*)

c-------------------------------------
c compute the Taylor deformation parameter
c         the  inclination                     
c-------------------------------------

      if(     Istep.eq.1.and.Iread.eq.0  ! special treatment of 
     +   .and.Ishape.eq.1) then          ! the spherical shape 
                                         ! commend out if desired
       rmax = 1.0
       rmin = 1.0
       zmax = 1.0
       dxy  (istep) =   0.0
c      thmax(istep) =  45.0
c      thmin(istep) = 135.0
       thmax(istep) =   0.25
       thmin(istep) =   0.75

      else

       call taylor 
     +
     +   (npts,nxy,jxy
     +   ,cx,cy,cz                ! deformation and
     +   ,rmax,rmin,zmax          ! inclination
     +   ,dxy(Istep)
     +   ,thmax(Istep),thmin(Istep)
     +   )

       call inclination
     +
     +    (nxy,jxy
     +    ,cx,cy            ! inclination by moment of
     +    ,thmax(Istep)    ! inertia tensor eigenvalues
     +    ,thmin(Istep)
     +    )

      end if

      write (6,*)
      write (4,104) Kstep,time(Kstep),rmax,rmin,Dxy(Kstep)
     +             ,thmax(Kstep),thmin(Kstep)
     +             ,zmax

      write (6,113)  time(Kstep),rmax,rmin,zmax
      write (6,115)   dxy(Kstep)
      write (6,116) thmax(Kstep)
      write (6,117) thmin(Kstep)
      write (6,*)

c---------------------------------
c prepare for the triply-periodic
c Green's function
c---------------------------------

      if(Iflow.eq.3) then

        if(a21.gt.a11h)  a21 = a21-a11
        if(a21.lt.a11hm) a21 = a21+a11

         call sgf_3d_3p_ewald
     +
     +     (a11,a12,a13
     +     ,a21,a22,a23
     +     ,a31,a32,a33
     +     ,b11,b12,b13
     +     ,b21,b22,b23
     +     ,b31,b32,b33
     +     ,ew,tau
     +     )

          call sgf_3d_3p_qqq
     +
     +      (b11,b12,b13
     +      ,b21,b22,b23
     +      ,b31,b32,b33
     +      ,Max2
     +      ,ew
     +      )

        if(Ivs.eq.1) then

         call sgf_3d_3p_vvv
     +
     +       (b11,b12,b13
     +       ,b21,b22,b23
     +       ,b31,b32,b33
     +       ,Max2
     +       ,ew
     +       )

        end if

       end if

c----------------------------------
c compute the velocity at the nodes
c----------------------------------
c
c     write (6,*) 
c     write (6,*) " Computing velocity at nodes"
c     write (6,*) 

      call caps_3d_vel
     +
     +    (npts
     +    ,nelm
     +    ,mint
     +    ,NGL
     +    ,Idfl
     +    ,Elst
     +    ,Elstb
     +    ,crvmr
     +    ,Iflow
     +    ,Isym_xy
     +    ,nter
     +    ,tol
     +    ,Istop,u
     +    )

      if(Istop.eq.1) Go to 99

c------------------------------------
c RHEOLOGY:
c
c In the case of infinite shear flow,
c compute the effective stress tensor
c------------------------------------

      if(Iflow.eq.1.or.Iflow.eq.3) then

      call rheology
     +
     +    (nelm
     +    ,npts
     +    ,mint
     +    ,vlm
     +    ,u
     +    ,sxy(istep)
     +    ,sd1(istep),sd2(istep)
     +    )

      if(abs(shrt).gt.0.00) then
        cf = 1.0/shrt              ! divide by the shear rate
      else
        cf = 1.0D0
      end if

      if(Iflow.eq.1) cf = cf/vlm

      sxy(Istep) = sxy(Istep)*cf
      sd1(Istep) = sd1(Istep)*cf
      sd2(Istep) = sd2(Istep)*cf

      write (6,*)
      write (6,118) sxy(Istep)
      write (6,119) sd1(Istep)
      write (6,120) sd2(Istep)

      write (6,*)
      write (6,*) "--------------------"
      write (6,*)

      end if

c------------------------------------
c compute the tank-treading frequency
c------------------------------------

      ttfto = 0.0D0
      ttftn = 0.0D0

      if(Isym_xy.eq.1) then

c---
c pick up position and velocity
c in the xy plane
c---

        Do i=1,nxy
           m = jxy(i)
           x_xy(i) = p(m,1)
           y_xy(i) = p(m,2)
          ux_xy(i) = u(m,1)
          uy_xy(i) = u(m,2)
        End Do

        call ttf
     +
     +     (nxy
     +     ,x_xy,y_xy
     +     ,ux_xy,uy_xy
     +     ,ttfto
     +     ,ttftn
     +     )

        write (6,121) ttfto
        write (6,122) ttftn

      end if

c--------------------
c print cross-section
c--------------------

      if(Isym_xy.eq.1) then

        write (6,*)
        write (6,*)  "position and velocity at time:",time(Istep)
        write (6,*)

c       Do i=1,nxy
c         m = jxy(i)
c         write (6,100) i,(p(m,j),j=1,2),(u(m,j),j=1,3)
c       End Do

        if(Iprint_xy.eq.Nprint_xy) then

          if(Iflow.eq.3) then
           write (3,100) nxy,time(Istep),a21
          else
           write (3,100) nxy,time(Istep)
          end if

          Do i=1,nxy
           m = jxy(i)
           write (3,100) i,(p(m,j),j=1,2)
          End Do

          Iprint_xy = 0

        end if

c---
      end if

c----------------
c print all nodes
c----------------

      if(Iprint_xyz.eq.Nprint_xyz) then

         if(Iflow.eq.3) then
          write (8,100) npts,time(Istep),a21
         else
          write (8,100) npts,time(Istep)
         end if

         Do i=1,npts
           write (8,108) i,(p(i,j),j=1,3)
         End Do

        ! print Matlab file (file 1)

        Do k=1,nelm
          call printel (k,Index,c)  ! print in file "caps_3d.net"
        End Do

        write (1,*) None

         Iprint_xyz = 0

      end if

c-----------------
c time integration
c-----------------

      if(Nstep.eq.0) Go to 99

c----------------------
      if(IRK.eq.2) then    ! second order

       Do i=1,npts
         ps(i,1) = p(i,1)   ! save position
         ps(i,2) = p(i,2)
         ps(i,3) = p(i,3)
       End Do

      end if
c----------------------

      Do i=1,npts

        if(move.eq.0) then      ! points move with total velocity

          Umove(i) = u(i,1)
          Vmove(i) = u(i,2)
          Wmove(i) = u(i,3)

        else if(move.eq.1) then  ! points move with normal velocity

          Uvel  = u(i,1)*vna(i,1)   ! projection of velocity
     +          + u(i,2)*vna(i,2)   ! onto the normal vector
     +          + u(i,3)*vna(i,3)

          Umove(i) = Uvel*vna(i,1)
          Vmove(i) = Uvel*vna(i,2)
          Wmove(i) = Uvel*vna(i,3)

        end if

c---
c flow in the presence of a wall
c---

        if(Iflow.eq.2) then            ! add tangential component
          if(move.eq.1) then           ! of centroid velocity

           Cntr_vel = shrt*(cy-wall)
           Umove(i) = Umove(i) + (1.0-vna(i,1)*vna(i,1))*Cntr_vel
           Vmove(i) = Vmove(i) + (   -vna(i,1)*vna(i,2))*Cntr_vel
           Wmove(i) = Wmove(i) + (   -vna(i,1)*vna(i,3))*Cntr_vel

          end if
        end if

c---
c advance in time
c---

        p(i,1) = p(i,1) + Dt*Umove(i)
        p(i,2) = p(i,2) + Dt*Vmove(i)
        p(i,3) = p(i,3) + Dt*Wmove(i)

c       write (6,100) i,(u(i,j),j=1,3) 

      End Do

c---
c Enforce symmetry if required
c---

      if(Isym_xy.eq.1) then
        Do i = 1,nxy
          p(jxy(i),3) = 0.0D0
        End Do
      end if

      if(Iflow.eq.3) a21 = a21 + Dt*shrt*a22    ! Triply-periodic flow

c------------
c end of RK1
c------------

       if(IRK.eq.1) Go to 87

c-----
c RK2: save the velocity
c------

      Do i=1,npts
        Usave(i) = Umove(i)
        Vsave(i) = Vmove(i)
        Wsave(i) = Wmove(i)
      End Do

c---
c Geometry at intermediate step
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
     +    (p(i1,1),p(i1,2),p(i1,3)
     +    ,p(i2,1),p(i2,2),p(i2,3)
     +    ,p(i3,1),p(i3,2),p(i3,3)
     +    ,p(i4,1),p(i4,2),p(i4,3)
     +    ,p(i5,1),p(i5,2),p(i5,3)
     +    ,p(i6,1),p(i6,2),p(i6,3)
     +    ,alpha(k),beta(k),gamma(k)
     +    )

      End Do

      call elm_geom 
     +
     +    (nelm,npts,mint
     +    ,area,vlm
     +    ,cx,cy,cz
     +    )

c----
c Prepare again for the
c triply-periodic Green's function
c-----

      if(Iflow.eq.3) then

        if(a21.gt.a11h)  a21 = a21-a11
        if(a21.lt.a11hm) a21 = a21+a11

         call sgf_3d_3p_ewald
     +
     +       (a11,a12,a13
     +       ,a21,a22,a23
     +       ,a31,a32,a33
     +       ,b11,b12,b13
     +       ,b21,b22,b23
     +       ,b31,b32,b33
     +       ,ew,tau
     +       )

          call sgf_3d_3p_qqq
     +
     +       (b11,b12,b13
     +       ,b21,b22,b23
     +       ,b31,b32,b33
     +       ,Max2
     +       ,ew
     +       )

        if(Ivs.eq.1) then

         call sgf_3d_3p_vvv
     +
     +       (b11,b12,b13
     +       ,b21,b22,b23
     +       ,b31,b32,b33
     +       ,Max2
     +       ,ew
     +       )
        end if

      end if

c---
c second velocity evaluation
c---

      call caps_3d_vel
     +
     +    (npts
     +    ,nelm
     +    ,mint
     +    ,NGL
     +    ,Idfl
     +    ,Elst
     +    ,Elstb
     +    ,crvmr
     +    ,Iflow
     +    ,Isym_xy
     +    ,nter
     +    ,tol
     +    ,Istop
     +    ,u
     +    )

      if(Istop.eq.1) Go to 99

c---
c second step in RK2
c
c see previous comments for an explanation 
c of the individual steps
c---

      Dth = 0.5D0*Dt

      Do i=1,npts

        if(move.eq.0) then

          Umove(i) = u(i,1)
          Vmove(i) = u(i,2)
          Wmove(i) = u(i,3)

        else if(move.eq.1) then

          Uvel  = u(i,1)*vna(i,1)
     +          + u(i,2)*vna(i,2)
     +          + u(i,3)*vna(i,3)

          Umove(i) = Uvel*vna(i,1)
          Vmove(i) = Uvel*vna(i,2)
          Wmove(i) = Uvel*vna(i,3)

        end if

        if(Iflow.eq.2) then

          if(move.eq.1) then 

           Cntr_vel = shrt*(cy-wall)

           Umove(i) = Umove(i) + (1.D0-vna(i,1)*vna(i,1))*Cntr_vel
           Vmove(i) = Vmove(i) + (    -vna(i,1)*vna(i,2))*Cntr_vel
           Wmove(i) = Wmove(i) + (    -vna(i,1)*vna(i,3))*Cntr_vel

          end if
        end if

        p(i,1) = ps(i,1) + Dth*(Usave(i)+Umove(i))
        p(i,2) = ps(i,2) + Dth*(Vsave(i)+Vmove(i))
        p(i,3) = ps(i,3) + Dth*(Wsave(i)+Wmove(i))

      End Do

c-----------
c end of RK2
c-----------

  87  Continue

c-------------------
c end of a time step
c-------------------

c------------------------------
c enforce symmetry if required
c------------------------------

      if(Isym_xy.eq.1) then
        Do i=1,nxy
         p(jxy(i),3) = 0.0D0
        End Do
      end if

c------------------------
c reset counters and time
c------------------------

      Kstep = Kstep+1
      Istep = Istep+1

      Iprint_xy  = Iprint_xy+1
      Iprint_xyz = Iprint_xyz+1

      time(Istep) = time(Istep-1) + Dt

      if(Istep.le.Nstep) Go to 97

      Go to 90     ! return for another step

c-.-.-.-.-.-.-.-.-.-.-
c simulation has ended
c-.-.-.-.-.-.-.-.-.-.-

  99  Continue

c-------------------
c record final shape
c-------------------

      open (11,file="caps_3d.rst")

      write (8,*) npts,time(Kstep)," restart data"
      write (11,*) npts,time(Kstep)," restart data"

      Do i=1,npts
        write (8,108) i,p(i,1),p(i,2),p(i,3)
        write (11,108) i,p(i,1),p(i,2),p(i,3)
      End Do

      write (8,100) Null
      write (11,100) Null

      close (11)

c-------------------
c record diagnostics
c-------------------

      write (3,*) Null
      write (4,*) Null

      write (4,*) Kstep," time, surf_ar, vol, cx, cy, cz"

      Do i=1,Kstep
        write (3,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (4,104) i,time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (6,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
        write (8,106)   time(i),ars(i),vol(i),crx(i),cry(i),crz(i)
      End Do

      write (8,*)
      write (3,*)
      write (4,*) Null
      write (4,*) Kstep," time, D, thmax, thmin"
      write (6,*)

      Do i=1,Kstep
        write (3,107)   time(i),dxy(i),thmax(i),thmin(i)
        write (4,104) i,time(i),dxy(i),thmax(i),thmin(i)
        write (6,107)   time(i),dxy(i),thmax(i),thmin(i)
        write (8,107)   time(i),dxy(i),thmax(i),thmin(i)
      End Do

      write (4,*) Null

c---
c record rheology
c---

      if(Iflow.eq.1) then

        write (3,*)
        write (4,*) Kstep," time, sxy,sd1, sd2"
        write (6,*)
        write (8,*)

        Do i=1,Kstep
          write (3,109)   time(i),sxy(i),sd1(i),sd2(i)
          write (4,104) i,time(i),sxy(i),sd1(i),sd2(i)
          write (6,109)   time(i),sxy(i),sd1(i),sd2(i)
          write (8,109)   time(i),sxy(i),sd1(i),sd2(i)
        End Do

        write (4,*) Null

      end if

c---
c record run parameters
c---

      write (3,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,req,boa,coa,cxs,cys,czs
     +             ,vs1,vs2,srft
     +             ,shrt
     +             ,Elst,Elstb
     +             ,crvmr

      write (4,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,req,boa,coa,cxs,cys,czs
     +             ,vs1,vs2,srft
     +             ,shrt
     +             ,Elst,Elstb
     +             ,crvmr

      write (6,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,req,boa,coa,cxs,cys,czs
     +             ,vs1,vs2,srft
     +             ,shrt
     +             ,Elst,Elstb
     +             ,crvmr

      write (8,205) ndiv,phi1,phi2,phi3
     +             ,Iflow,wall
     +             ,req,boa,coa,cxs,cys,czs
     +             ,vs1,vs2,srft
     +             ,shrt
     +             ,Elst,Elstb
     +             ,crvmr

      write (3,208) mint,NGL,dt,irk,move,norm,isym_xy
      write (4,208) mint,NGL,dt,irk,move,norm,isym_xy
      write (6,208) mint,NGL,dt,irk,move,norm,isym_xy
      write (8,208) mint,NGL,dt,irk,move,norm,isym_xy

c--------------------------
c print Matlab file (file 1)
c--------------------------

        Do k=1,nelm
          call printel (k,Index,c)  ! print in file "caps_3d.net"
        End Do

        write (1,*) Null

c---
c wrap up
c---

      close (1)
      close (2)
      close (3)
      close (4)
      close (8)

c-----
c Done
c-----

  100 Format (1x,i4,10(1x,f12.5))
  101 Format (10(1x,f12.5))
  102 Format (10(1x,f10.6))
  103 Format (1x,i4,10(1x,f10.5))
  104 Format (1x,i4,10(1x,f8.5))

  105 Format (" Step ",i4," out of ",i4,"; time :",F15.10)
  106 Format (' T=',F7.3,' S=',F10.7,' V=',F10.7
     +       ,' X=',F8.5,' Y=',F 8.5,' Z=',F8.5
     +       )
  107 Format (' T=',F7.3,' D=',F10.7,' thmax=',F10.4,
     +                               ' thmin=',F10.4)
  108 Format (1x,i4,100(1x,f15.10))
  109 Format (' T=',F7.3,' sxy=',F10.6,' sd1=',F10.6,
     +                                 ' sd2=',F10.6)

  110 Format (" Surface Area :",F15.10)
  111 Format (" Volume       :",F15.10)
  112 Format (" Centroid     :",3(F15.10))
  113 Format ("Time=",F12.5," Axes:",3(1x,F12.5))

  115 Format (" Deformation  :",F15.10)
  116 Format (" Max Incl     :",F10.4)
  117 Format (" Min Incl     :",F10.4)
  118 Format (" Eff shear st :",F10.4)
  119 Format (" Eff first nsd:",F10.4)
  120 Format (" Eff sec   nsd:",F10.4)

  121 Format (" tank-Treading frequency based on",
     + " total velocity     :",F15.10)

  122 Format (" tank-Treading frequency based on",
     + " tangential velocity:",F15.10)

  200 Format(100(1x,f5.3))

  205 Format (/
     +       ,' ndiv   = ',I2,/
     +       ,' phi1   = ',F7.4,/
     +       ,' phi2   = ',F7.4,/
     +       ,' phi3   = ',F7.4,/,/
     +       ,' Iflow  = ',I2,/
     +       ,' wall   = ',F7.4,/,/
     +       ,' eq rad = ',F7.4,/
     +       ,' b/a    = ',F7.4,/
     +       ,' c/a    = ',F7.4,/
     +       ,' cx     = ',F7.4,/
     +       ,' cy     = ',F7.4,/
     +       ,' cz     = ',F7.4,/,/
     +       ,' vs1    = ',F10.5,/
     +       ,' vs2    = ',F10.5,/
     +       ,' srft   = ',F7.4,/,/
     +       ,' shrt   = ',F7.4,/
     +       ,' Elst   = ',F12.6,/
     +       ,' Elstb  = ',F12.6,/
     +       ,' crvmr  = ',F12.6,/
     +       )

 208  Format
     +       (' mint   = ',I1,/
     +       ,' NGL     = ',I1,/
     +       ,' DT     = ',F8.6,/
     +       ,' RUNGE  = ',I1,/
     +       ,' MOVE   = ',I1,/
     +       ,' Norm   = ',I1,/
     +       ,' Isym_xy= ',I1,/
     +       )

      stop
      end
