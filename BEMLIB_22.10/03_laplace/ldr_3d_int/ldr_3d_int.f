      program ldr_3d_int

c=========================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------
c Solution of the Dirichlet problem for
c Laplace's equation in the interior of
c a closed three-dimensional region
c
c This program solves the equation:
c
c     Laplacian(f) = 0.0
c
c subject to Dirichlet boundary conditions
c for the desired function f
c at the external (closed) boundary
c
c The solution is represented by a dipole
c distribution over the exterior boundary
c
c SYMBOLS:
c -------
c
c f(i):      boundary value of f at the ith node
c q(i):      normal derivative of f at the ith node
c            q = n.grad(f)
c            where n is the outward normal vector
c dpl(i):    density of the dipole distribution at the ith node
c dlp_pv(i): principal value of the double-layer potential
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension      p(1026,3)
      Dimension     ne(1026,7)
      Dimension    vna(1026,3)
      Dimension      f(1026),dpl(1026)
      Dimension dlp_pv(1026)
      Dimension sgrad_f(1026,3),gradn_f(1026,3),gradn_f_n(1026)
      Dimension    zeta(1026,3),slp3(1026,3)

      Dimension       n(512,6), nbe(512,3)
      Dimension   alpha(512),  beta(512),gamma(512)
      Dimension    arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension  crvmel(512)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/var1/Iflow,Ign,wall

      common/geo1/arel
      common/geo2/vna
      common/geo6/crvmel

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      piq = 0.25D0*pi
      pih = 0.50D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi
      pi6 = 6.0D0*pi
      pi8 = 8.0D0*pi

      null = 0
      Nfour  = 4
      Nseven = 7
      oot  = 1.0D0/3.0D0

c-----------
c input data
c-----------

      write (6,*)
      write (6,*) ' Enter :'
      write (6,*)
      write (6,*) ' 1 to read them from file: ldr_3d_int.dat'
      write (6,*) ' 2 to type the data into the keyboard'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*)   Ienrd

c------------------------
      if(Ienrd.eq.1) then
c------------------------

       open (2,file="ldr_3d_int.dat")

       read (2,*) Ioctaicos
       read (2,*) ndiv
       read (2,*)
       read (2,*) boa,coa
       read (2,*) req
       read (2,*) cxp,cyp,czp
       read (2,*)
       read (2,*) mint
       read (2,*) NGL
       read (2,*)
       read (2,*) Nter
       read (2,*) tol
       read (2,*)
       read (2,*) Iflow
       read (2,*)
       read (2,*) wall
       read (2,*) Ign

      close (2)

c---------
      else
c---------

  91  Continue

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to discretize an octahedron"
      write (6,*) " 2 to discretize an icosahedron"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ioctaicos

      if(Ioctaicos.eq.0) Go to 99

      write (6,*)
      write (6,*) " Select the level of triangulation"
      write (6,*) " Choose from 0, 1, 2, 3; 99 to quit"
      write (6,*) " ----------------------------------"
      read  (5,*) ndiv

      if(ndiv.eq.99) Go to 99


      if((ndiv.lt.0).or.(lt.gt.3))  then
        write (6,*) 'out of range; please try again'
        Go to 91
      end if

      write (6,*)
      write (6,*) " The x,y,z semi-axes of the ellipsoid are: a,b,c"
      write (6,*)
      write (6,*) " Enter the axes ratios b/a and c/a"
      write (6,*) " ---------------------------------"
      read  (5,*) boa,coa

      write (6,*)
      write (6,*) " Enter the ellipsoid equivalent radius"
      write (6,*) " -------------------------------------"
      read  (5,*) req


      write (6,*)
      write (6,*) " Enter the coordinates of the"
      write (6,*) " center of the ellipsoid"
      write (6,*) " -----------------------"
      read  (5,*) cxp,cyp,czp

      write (6,*)
      write (6,*) " Will use the m-point gauss-triangle rule"
      write (6,*) " for integration over each triangle"
      write (6,*)
      write (6,*) " Please enter m"
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*) "----------------"
      read  (5,*) mint

      write (6,*)
      write (6,*) " Enter the order of the quadrature"
      write (6,*) " for the singular Gauss-Legendre"
      write (6,*) " integration"
      write (6,*) " -----------"
      read  (5,*) NGL

      write (6,*)
      write (6,*) " Enter the maximum number of interations"
      write (6,*) " in solving the deflated integral equation"
      write (6,*) "------------------------------------------"
      read  (5,*) Nter

      write (6,*)
      write (6,*) " Enter the tolerance"
      write (6,*) " -------------------"
      read  (5,*) tol

      write (6,*)
      write (6,*) " Choose the Green's function"
      write (6,*)
      write (6,*) " Enter 1 for free space"
      write (6,*) "       2 for domain bounded by a plane wall"
      write (6,*) " ------------------------------------------"
      read  (5,*) Iflow

      if(Iflow.eq.2) then

        write (6,*) " The wall is located at x=wall"
        write (6,*) " Please enter: wall"
        write (6,*) " ------------------"
        read  (5,*) wall

        write (6,*) " Enter 1 for the Green function"
        write (6,*) "       2 for the Neumann function"
        write (6,*) " --------------------------------"
        read  (5,*) Ign

      end if

c-----------
      end if
c-----------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for f = 1.0'
      write (6,*) ' 2 for f = x'
      write (6,*) ' 3 for f = x**2'
      write (6,*) ' 4 for f = x*y*z'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) Ibc

      if(Ibc.eq.0) Go to 99

c-----------------
c open output file
c-----------------

      open (4,file="ldr_3d_int.out")

c----------------------------
c triangulate the unit sphere
c----------------------------

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

c---
c Expand to specified radius;
c Deform to an ellipsoid;
c Translate center to specified position.
c---

      scale = req/(boa*coa)**oot

      Do i=1,npts
        p(i,1) = scale*    p(i,1) + cxp
        p(i,2) = scale*boa*p(i,2) + cyp
        p(i,3) = scale*coa*p(i,3) + czp
      End Do

c---------------------------------
c specify the boundary values of f
c---------------------------------


      Do i=1,npts

       if(Ibc.eq.1) then
           f(i) = 1.0
       else if(Ibc.eq.2) then
           f(i) = p(i,1)
       else if(Ibc.eq.3) then
           f(i) = p(i,1)**2
       else if(Ibc.eq.4) then
           f(i) = p(i,1)*p(i,2)*p(i,3)
       end if

      End Do

c-----------------------------
c read the triangle quadrature
c-----------------------------

      call gauss_trgl
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

c----------------------------
c compute alpha, beta, gamma,
c for the xi-eta mapping
c over each triangle
c---------------------------

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

c------------------------------------------------
c compute:
c
c  arel:            surface area of the individual elements
c  xmom, ymom,zmom: x, y, and z surface moments
c                   over each element
c  area:            total surface area
c  vlm:             volume
c  crvmel:           mean curvature over each element
c  vna:             averaged normal vector at the nodes
c------------------------------------------------

      call elm_geom
     +
     +   (nelm,npts
     +   ,mint
     +   ,xmom,ymom,zmom
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c---
c normalize surface area and volume
c---

      arean = area/(pi4*req**2)         ! normalize
      vlmn  = 3.0D0*vlm/(pi4*req**3)    ! normalize

      write (6,*)
      write (6,110) arean
      write (6,111) vlmn
      write (6,*)
      write (6,800) cx,cy,cz

c----------------------
c initialize the dipole (typical)
c----------------------

      Do i=1,npts
       dpl(i) = 0.9876
c      dpl(i) = 0.1
      End Do

c--------------------------------------------
c Iterative solution of the integral equation
c for the dipole strength
c--------------------------------------------


      Iter = 0

 71   Continue  ! target for iterations

      Iter = Iter+1

c-----------------------------------
c compute the double-layer potential
c-----------------------------------

      call ldlp_3d
     +
     +  (npts
     +  ,nelm
     +  ,mint
     +  ,dpl
     +  ,dlp_pv
     +  )

c     Do i=1,npts
c      write (6,100) i,dpl(i),dlp_pv(i)
c     End Do
c     pause

c-----------------------------
c compute the surface integral
c of the dipole density
c-----------------------------

      call srf_int_3d
     +
     +   (nelm,npts
     +   ,mint
     +   ,dpl
     +   ,dpl_sint
     +   )

c-----------------------------
c update and compute the error
c-----------------------------

      error = 0.0D0
     
      Do i=1,npts

       save = dpl(i)

       dpl(i) = 2.0D0*(dlp_pv(i)-f(i)) + dpl_sint/area

       error = error + (save-dpl(i))**2

c      write (6,100) i,f(i),dpl(i),dlp_pv(i)

      End Do

      error = Dsqrt(error)/npts

      write (6,456) Iter,dpl_sint,error

      if(error.lt.tol) Go to 711
      if(Iter.lt.Nter) Go to 71

      write (6,*)
      write (6,*) " One more iteration batch ?"
      write (6,*)
      write (6,*) " 1 for Yes, 0 for NO"
      write (6,*) "--------------------"
      read  (5,*) more

      if(more.eq.1) then
       Iter = 0
       Go to 71
      end if

c-----------------------------------------------------
c recover the solution
c using equation (10.7.20) of Pozrikidis (1997, p.514)
c-----------------------------------------------------

 711   Continue

      write (6,*)
      write (6,*) "converged solution: x-position, dipole"
      write (6,*)

      Do i=1,npts

       dpl(i) = dpl(i)-0.5D0*dpl_sint/area

       write (4,100) i,p(i,1),dpl(i)
       write (6,100) i,p(i,1),dpl(i)

      End Do

c---------------------
c generate a Matlab file
c---------------------

      open (1,file="ldr_3d_int.net")

      index = 1  ! 6-node triangles
      index = 2  ! 3-node triangles

      if(Index.eq.1) then          ! 6-node triangles
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*)   nelm
      else if(Index.eq.2) then     ! 3-node triangles
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*)  4*nelm
      end if

      write (1,*) Iflow
      write (1,*) wall

      Do k=1,nelm
       call printel (k,Index,dpl)  ! print in file "caps_3d.net"
      End Do

      close (1)

c----------------
c post-processing
c----------------

 97   Continue

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to evaluate the solution at a point"
      write (6,*) " 2 to compute the surface gradient"
      write (6,*) " 3 to compute the derivative normal to"
      write (6,*) "   the boundary and its surface integral"
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c==========================
      if(menu.eq.1) then
c==========================

      write (6,*)
      write (6,*) " Enter the coordinates of the evaluation point"
      write (6,*) "----------------------------------------------"
      read  (5,*) xp,yp,zp

      call ldlpp_3d
     +
     +  (npts,nelm
     +  ,mint
     +  ,dpl
     +  ,xp,yp,zp
     +  ,ptl
     +  )

      write (6,*)
      write (6,*) "The solution is:"
      write (6,101) ptl

c===========
      end if
c===========

c=======================
      if(menu.eq.2) then   ! compute the surface gradient
c=======================

      call sgrad_3d
     +
     +  (npts,nelm
     +  ,f
     +  ,sgrad_f
     +  )

      write (6,*)
      write (6,*) "x-position - surface gradient"
      write (6,*)

      Do i=1,npts
       write (6,100) i,p(i,1),sgrad_f(i,1),sgrad_f(i,2),sgrad_f(i,3)
      End Do

c============
      end if
c============

c=======================
      if(menu.eq.3) then  ! normal derivative and its integral
c=======================

c---
c compute the strength of the
c vortex sheet
c from the dipole strength
c using the equation:
c
c zeta = n x grad(dpl)
c---

      call vs_3d_circ_vort
     +
     +   (npts,nelm
     +   ,dpl
     +   ,zeta
     +   )

c---
c call the Gauss--Legendre quadrature
c---

      call gauss_leg (NGL,zz,ww)

c---
c compute the vector potential
c at the nodes as a single-layer
c potential
c---

      call lslp3_3d
     +
     +   (npts,nelm
     +   ,mint,NGL
     +   ,zeta
     +   ,slp3
     +   )

c---
c compute the normal component
c of the curl of slp3
c
c gradn_f = n (n.gradf)
c---

      call vs_3d_curl
     +
     +    (npts,nelm
     +    ,slp3
     +    ,gradn_f
     +    )

c---
c compute the projection
c of the curl of slp3
c onto the normal vector
c
c gradn_f_n = n.gradn_f
c---

      write (6,*)
      write (6,*) "x-position - magnitude of n.grad(f)"
      write (6,*)

      Do i=1,npts

       gradn_f_n(i) = gradn_f(i,1)*vna(i,1)
     +               +gradn_f(i,2)*vna(i,2)
     +               +gradn_f(i,3)*vna(i,3)

       write (6,100) i,p(i,1),gradn_f_n(i)

      End Do

c---
c compute the surface integral
c of gradn_f_n
c---

      call srf_int_3d
     +
     +    (nelm,npts,mint
     +    ,gradn_f_n
     +    ,flow
     +    )

      write (6,*)
      write (6,457) flow

c--
c record for visualization
c--

      open (1,file="ldr_3d_int_1.net")

      Index = 1  ! 6-node triangles
      Index = 2  ! 3-node triangles

      if(Index.eq.1) then          ! 6-node triangles
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*)   nelm
      else if(Index.eq.2) then     ! 3-node triangles
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*)  4*nelm
      end if

      write (1,*) Iflow
      write (1,*) wall

      Do k=1,nelm
       call printel (k,Index,gradn_f_n)  ! print in file "caps_3d.net"
      End Do

      close (1)

c===========
      end if
c===========

      Go to 97    ! return to the Menu

c-----
c done
c-----

 99   Continue

      close (4)

  100 Format (1x,i3,10(1x,f10.5))
  101 Format (f15.10)

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)

  456 Format ("Iteration: ",I3," Projection:",f15.8," Error: ",f20.15)
  457 Format ("Flow across the boundary:",f15.10)

  800 Format(" centroid:",3(1x,f10.5))

      Stop
      End
