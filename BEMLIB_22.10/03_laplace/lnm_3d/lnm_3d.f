      program lnm_3d 

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c------------------------------------------
c Solution of Laplace's equation with the Neumann
c boundary condition in the interior or exterior
c of a closed three-dimensional surface
c
c This program solves the equation:
c 
c Laplacian(f) = 0.0
c
c subject to the Neumann boundary condition
c for the normal derivative
c dfdn = n.grad(f) over the closed surface
c
c where n is the outward normal vector
c
c The solution is based on the standard boundary-integral
c representation (Green's third identitity)
c
c SYMBOLS:
c -------
c
c f(i):   specified boundary value of f at the ith node
c
c dfdn(i):   normal derivative of f at the ith node
c            where n is the outward normal vector
c            dfdn = n.grad(f)
c
c dlp(i): double-layer potential at the ith node
c
c q(i): auxiliary variable used for impulses
c       of the double-layer potential
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)
      Dimension   f(1026),q(1026)
      Dimension dlp(1026)

      Dimension     n(512,6), nbe(512,3)
      Dimension alpha(512),  beta(512),gamma(512)
      Dimension  arel(512),  xmom(512),  ymom(512),zmom(512)

      Dimension amat(1026,1026),rhs(1026),dfdn(1026)

      Dimension zz(20),ww(20)
      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe

      common/albega/alpha,beta,gamma

      common/geo1/arel
      common/geo2/vna

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      piq = 0.25D0 *pi
      pih = 0.50D0 *pi
      pi2 = 2.00D0 *pi
      pi4 = 4.00D0 *pi
      pi6 = 6.00D0 *pi
      pi8 = 8.00D0 *pi

      Null = 0
      Nfour = 4
      Nseven = 7
      oot  = 1.0D0/3.0D0

c-----------
c Input data
c-----------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to read them from file: lnm_3d.dat'
      write (6,*) ' 2 to type the data into the keyboard'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*)   Ienrd

      if(Ienrd.eq.0)  Go to 99

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (2,file="lnm_3d.dat")

       read (2,*) int_ext
       read (2,*)
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
       read (2,*) Iflow

      close (2)

c---------
      else
c---------

  91  Continue

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for the interior problem"
      write (6,*) " 2 for the exterior problem"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      write (6,*)
      read  (5,*) int_ext

      if(int_ext.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to discretize an octahedron"
      write (6,*) " 2 to discretize an icosahedron"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Ioctaicos

      if(Ioctaicos.eq.0) Go to 99

      write (6,*) " Select the level of triangulation"
      write (6,*) " Enter 0, 1, 2, 3; 99 to quit"
      write (6,*) " ----------------------------"
      read  (5,*) ndiv

      if(ndiv.eq.99) Go to 99

      if(ndiv.gt.3) then
        write (6,*) 'out of range; please try again.'
        Go to 91
      End If

      write (6,*)
      write (6,*) " Enter the axes ratios b/a and c/a"
      write (6,*) " ---------------------------------"
      read  (5,*) boa,coa

      write (6,*)
      write (6,*) " Enter the equivalent radius"
      write (6,*) " ---------------------------"
      read  (5,*) req

      write (6,*)
      write (6,*) " Enter the coordinates of the "
      write (6,*) " center of the ellipsoid"
      write (6,*) " ------------------------------"
      read  (5,*) cxp,cyp,czp

      write (6,*)
      write (6,*) " Will use the m-point gauss-triangle rule"
      write (6,*) " for integrating over each triangle"
      write (6,*)
      write (6,*) " Please enter m"
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*) "----------------"
      read  (5,*) mint

      write (6,*)
      write (6,*) " Enter the number of Gauss-Legendre"
      write (6,*) " nodes for integrating the single-layer"
      write (6,*) " potential in polar coordinates"
      write (6,*) " ------------------------------"
      read  (5,*) NGL

      write (6,*)
      write (6,*) " Choose the flow:"
      write (6,*)
      write (6,*) " Enter 1 for the free-space"
      write (6,*) " ---------------------------"
      read  (5,*) Iflow

c-----------
      end if
c-----------

c--------
c Warning
c--------

      if(int_ext.eq.1) then
        write (6,*)
        write (6,*) 'WARNING'
        write (6,*)
        write (6,*) 'You are attempting the solve the interior'
        write (6,*) 'Neumann problem; the exact solution either'
        write (6,*) 'does not exist, or it is not unique'
        write (6,*) 
        write (6,*) 'Proceed at your own risk'
        write (6,*) 
        pause
      end if

c--------
c boundary data
c--------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for dfdn = 1.0'
      write (6,*) ' 2 for dfdn = x'
      write (6,*) ' 3 for dfdn = x**2'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) Ibc

      if(Ibc.eq.0) Go to 99

c-----
c trap
c-----

      if(Iflow.ne.1) then
        write (6,*) "lnm_3d: only the free-space green's function"
        write (6,*) "        is implemented; our apologies"
        stop
      end if

c------------------
c open output files
c------------------

      open (4,file="lnm_3d.out")

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
      write (6,*) " lnm_3d: number of points:   ",npts
      write (6,*)
      write (6,*) " lnm_3d: number of elements: ",nelm
      write (6,*)

c---
c expand to specified radius;
c deform to an ellipsoid;
c translate center
c to specified position;
c---

      scale = req/(boa*coa)**oot

      Do i=1,npts
        p(i,1) = scale*    p(i,1) + cxp
        p(i,2) = scale*boa*p(i,2) + cyp
        p(i,3) = scale*coa*p(i,3) + czp
      End Do

c------------------------------------
c specify the boundary values of dfdn
c------------------------------------

      Do i=1,npts

       if(Ibc.eq.1) then
         dfdn(i) = 1.0D0
       else if(Ibc.eq.2) then
          dfdn(i) = p(i,1)
       else
          dfdn(i) = p(i,1)**2
       end if

      End Do

c---------------------
c read the quadratures
c---------------------

      call gauss_leg (NGL,zz,ww)

      call gauss_trgl
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

c---------------------------------------------
c compute the coefficients alpha, beta, gamma,
c for the quadratic xi-eta mapping
c of each element
c---------------------------------------------

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
c  vlm:  surface area of the individual elements
c  xmom,yomo,zmom:  x, y, and z surface moments
c                      over each element
c  area:   total surface area and volume
c  crvmel: mean curvature over each element
c  vna:    average normal vector at the nodes
c------------------------------------------------

      call elm_geom
     +
     +   (nelm,npts,mint
     +   ,xmom,ymom,zmom
     +   ,area,vlm
     +   ,cx,cy,cz
     +   )

c---
c normalize surface area and volume
c---

      area = area/(pi4*req**2)          ! normalize
      vlm  = 3.0D0*vlm /(pi4*req**3)    ! normalize

      write (6,*)
      write (6,110) area
      write (6,111) vlm
      write (6,*)

      write (6,800) cx,cy,cz

c--------------------------
c assemble the linear system
c--------------------------

c----------------------------
c compute the right-hand side
c of the linear system
c defined in terms of the
c single-layer potential
c----------------------------

      write (6,*)
      write (6,*)  " lnm_3d: computing the right-hand side "
      write (6,*)  "         of the linear system"
      write (6,*)

       call lslp_3d
     +
     +    (npts
     +    ,nelm
     +    ,mint,NGL
     +    ,dfdn
     +    ,Iflow
     +    ,rhs
     +    )

      write (6,*) 'right-hand side computed'
      write (6,*)

c---
c display the rhs
c---

      write (6,*)
      write (6,*)  " lnm_3d: right-hand side of the linear system:"
      write (6,*)

      Do i=1,npts
       write (6,100) i,rhs(i)
      End Do

c-----------------------------
c compute the influence matrix
c by the method of impulses
c-----------------------------

c---
c initialize
c---

      Do i=1,npts
       q(i) = 0.0D0
      End Do

c---
c proceed with the impulses
c---

      Do j=1,npts

       write (6,*) "Computing column ",j, " of the influence matrix"

       q(j) = 1.0D0   ! impulse

       call ldlp_3d
     +
     +   (npts
     +   ,nelm
     +   ,mint
     +   ,q
     +   ,Iflow
     +   ,dlp
     +   )

       Do k=1,npts
        amat(k,j) = dlp(k)
       End Do

       q(j) = 0.0D0    ! reset

      End Do

c---
c complete the influence matrix
c---

      if(int_ext.eq.1) then
         cf = 0.5D0                ! interior
      else if(int_ext.eq.2) then
         cf =-0.5D0                ! exterior
      end if

      Do i=1,npts
       amat(i,i) = amat(i,i) + cf
      End Do

c---
c testing
c---

c     Do j=1,npts
c      test = 0.0
c      Do k=1,npts
c       test = test+amat(j,k)
c      End Do
c      write (6,*) j,test
c     End Do

c--------------------------

      write (6,*) 'linear system assembled'
      write (6,*)

c----------------------------------
c done generating the linear system
c solve by gauss elinimation
c----------------------------------

      Isym_gel = 0   ! system is not symmetric
      Iwlpvt   = 1   ! pivoting enabled

      call gel
     +
     +   (npts,amat
     +   ,rhs
     +   ,f
     +   ,Isym_gel
     +   ,Iwlpvt
c    +   ,l,u
c    +   ,det
     +   ,Istop
     +   )

c--------------------------------
c record and display the solution
c--------------------------------

      write (6,*)
      write (6,*) "node coordinates and solution"
      write (6,*)

      write (4,100) npts

      Do i=1,npts
       write (4,100) i,p(i,1),p(i,2),p(i,3),f(i)
       write (6,100) i,p(i,1),p(i,2),p(i,3),f(i)
      End Do

c--------------------
c print the triangles
c and visualize the solution
c for the normal derivative
c--------------------

      open (1,file="ldr_3d.net") ! output for matlab visualization

      Index = 1  ! 6-node triangles
      Index = 2  ! 3-node triangles

      if(Index.eq.1) then
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*) nelm
      else if(Index.eq.2) then
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*) 4*nelm
      end if

      write (1,*) Iflow

      Do k=1,nelm
        call printel (k,index,f)
      End Do

      close(1)

c----------------
c post processing
c----------------

 97   Continue

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to evaluate the solution at an arbitrary point"
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c-----------------------
      if(menu.eq.1) then
c-----------------------

      write (6,*)
      write (6,*) " Enter the coordinates of the point"
      write (6,*) " where the solution is to be evaluated"
      write (6,*) "--------------------------------------"
      read  (5,*) xp,yp,zp

      call lsdlpp_3d
     +
     +    (nelm
     +    ,mint
     +    ,f
     +    ,dfdn
     +    ,xp,yp,zp
     +    ,Iflow
     +    ,fp
     +    )

      if(int_ext.eq.1) fp = -fp

      write (6,*) "lnm_3d: solution:"

      write (6,101) fp

c-----------
      end if
c-----------

      Go to 97   ! return to menu

c-----
c done
c-----

 99   Continue

      write (4,*) Null
      close (4)

  100 Format (1x,i3,10(1x,f10.5))
  101 Format (f15.10)

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)

  800 Format(" centroid:",3(1x,f10.5))

      stop
      end
