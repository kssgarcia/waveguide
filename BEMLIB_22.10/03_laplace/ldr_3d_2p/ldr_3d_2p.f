      program ldr_3d_2p

c=============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c=============================================

c------------------------------------------
c Solution of Laplace's equation
c in a semi-infinite domain
c bounded by a doubly-periodic surface
c that is conformal to the xy plane
c with Dirichlet condition
c
c This program solves the equation:
c
c Laplacian(f) = 0.0,
c
c subject to a boundary condition for the function f
c at the doubly-periodic boundary
c
c The solution is represented in terms of a 
c doublye periodic harmonic double-layer potential.
c
c SYMBOLS:
c -------
c
c g(i):   boundary value of f at the ith node
c q(i):   normal derivative of f at the ith node
c         q = n.grad(f)
c         where n is the normal vector pointing upward
c
c dpl(i):  density of the dipole distribution at the ith node
c
c dlp_pv(i): principal value of the double-layer potential
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension       p(1090,3)
      Dimension      ne(1090,7)
      Dimension     vna(1090,3)
      Dimension   Iedge(1090,4)
      Dimension       g(1090),dpl(1090)
      Dimension    auxs(1090)
      Dimension    zeta(1090,3)
      Dimension    slp3(1090,3)
      Dimension  dlp_pv(1090)
      Dimension sgrad_f(1090,3),gradn_f(1090,3),gradn_f_n(1090)

      Dimension      n(512,6), nbe(512,3)
      Dimension  alpha(512),  beta(512),gamma(512)
      Dimension   arel(512),  xmom(512),  ymom(512),zmom(512)
      Dimension crvmel(512)

      Dimension zz(20),ww(20)
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
      common/geo6/crvmel

      common/zwl/zz,ww
      common/trq/xiq,etq,wq

      common/lgfi/method_lgf,Max1,Max2,Max3
      common/lgfr/a11,a12,a21,a22,b11,b12,b21,b22,ew,area_lgf

      common/pii/pi,piq,pih,pi2,pi4,pi6,pi8

c----------
c constants
c----------

      pi = 3.14159 265358 D0

      piq = 0.25D0*pi
      pih = 0.50D0*pi
      pi2 = 2.00D0*pi
      pi4 = 4.00D0*pi
      pi6 = 6.00D0*pi
      pi8 = 8.00D0*pi

      Null  = 0
      Nfour  = 4
      Nseven = 7

c---
c Input data
c---

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to read them from file: ldr_3d_2p.dat'
      write (6,*) ' 2 to type the data into the keyboard'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*)  Ienrd

      Ienrd = 1

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (2,file="ldr_3d_2p.dat")

       read (2,*) cud    ! coefficnet up/down
       read (2,*) grad_inf
       read (2,*)
       read (2,*) ndiv
       read (2,*)
       read (2,*) a11,a12
       read (2,*) a21,a22
       read (2,*)
       read (2,*) amp1,amp2
       read (2,*)
       read (2,*) mint
       read (2,*) NGL
       read (2,*)
       read (2,*) Method_lgf
       read (2,*) Max1,Max2,Max3
       read (2,*)
       read (2,*) Nter
       read (2,*) error_max

      close (2)

c---------
      else
c---------

  91  Continue

      write (6,*)
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) " 1 to select to upper semi-space"
      write (6,*) "-1 to select to lower semi-space"
      write (6,*) ' 0 to quit'
      write (6,*) " ---------"
      read  (5,*) cud

      if(abs(cud).le.0.000000001) Go to 99

      write (6,*)
      write (6,*) " Boundary condition at infinity is: "
      write (6,*) "    f -> grad_inf * z"
      write (6,*)
      write (6,*) " Enter the coefficient: grad_inf"
      write (6,*) " --------------------------------"
      read  (5,*) grad_inf

      write (6,*) " Select the level of triangulation"
      write (6,*)
      write (6,*) " Enter 0, 1, 2, 3; 99 to quit"
      write (6,*) " ----------------------------"
      read  (5,*) ndiv

      if(ndiv.eq.99) Go to 99

      if(ndiv.gt.3)  then
        write (6,*) "out of range; please try again"
        Go to 91
      end if

      write (6,*)
      write (6,*) " Enter the components of the base vectors"
      write (6,*) " in the xy plane"
      write (6,*) " ---------------"
      read  (5,*) a11,a12
      read  (5,*) a21,a22

      write (6,*)
      write (6,*) " Enter the amplitudes of the wall in the"
      write (6,*) " directions of the base vectors" 
      write (6,*) " ---------------------------------------"
      read  (5,*) amp1,amp2

      write (6,*)
      write (6,*) " Will use the m-point gauss-triangle rule"
      write (6,*) "      for integration over each triangle"
      write (6,*)
      write (6,*) " Please enter m"
      write (6,*)
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
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
      write (6,*) " Choose the method of computing the"
      write (6,*) "        doubly-periodic Laplace Green funcion"
      write (6,*)
      write (6,*) " Enter 1 for Fourier series"
      write (6,*) "       2 for Ewald sums"
      write (6,*) " ----------------------"
      read  (5,*) Method_lgf

      if(Method_lgf.eq.1) then

       write (6,*)
       write (6,*) " Enter the summation truncation limit: Max"
       write (6,*) " -----------------------------------------"
       read  (5,*) Max2

      else

       write (6,*)
       write (6,*) " Enter the summation truncation limits:"
       write (6,*) " Max1, Max2, Max3"
       write (6,*) " ----------------"
       read  (5,*) Max1,Max2,Max3

      end if

      write (6,*)
      write (6,*) "Enter the maximum number of interations"
      write (6,*) "in solving the deflated integral equation"
      write (6,*) "------------------------------------------"
      read  (5,*) Nter

      write (6,*)
      write (6,*) "Enter the maximum error (tolerance)"
      write (6,*) "------------------------------------"
      read  (5,*) error_max

c-----------
      end if
c-----------

      write (6,*)
      write (6,*) ' Choose the boundary data'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for g = 0.0'
      write (6,*) ' 2 for g = z'
      write (6,*) ' 3 for g = z**2'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) Ibc

      if(Ibc.eq.0) Go to 99

c-----------------
c open output file
c-----------------

      open (4,file="ldr_3d_2p.out")
      open (3,file="ldr_3d_2p.grid")

c----------------------------
c triangulate the unit square
c----------------------------

      call trgl6_sqr
     +
     +  (ndiv
     +  ,npts,nelm
     +  )

      write (6,*)
      write (6,*) "Number of points:   ",npts
      write (6,*)
      write (6,*) "Number of elements: ",nelm
      write (6,*)

c---
c Display the edge nodes
c         and their images
c---

      write (6,*)
      write (6,*) " Edge nodes and images"
      write (6,*)

      Do i=1,npts
       if(Iedge(i,1).ne.0) then
        write (6,139) i,(Iedge(i,j),j=1,Iedge(i,1)+1)
       End if
      End Do

c------------------------------------
c Deform the unit square to specified
c shape determined by the base vectors
c a_1 and a_2
c
c Impose a periodic perturbation 
c------------------------------------

      cell_area = a11*a22-a21*a12

      cf  = pi2/cell_area

      b11 =  cf*a22
      b12 = -cf*a21
      b21 = -cf*a12
      b22 =  cf*a11

c     write (6,100) npts,time(Istep)

      Do i=1,npts
        tmpx   = p(i,1)
        tmpy   = p(i,2)
        p(i,1) = a11*tmpx+a21*tmpy
        p(i,2) = a12*tmpx+a22*tmpy
        arg1   = b11*p(i,1)+b12*p(i,2)
        arg2   = b21*p(i,1)+b22*p(i,2)
        p(i,3) = amp1*cos(arg1)+amp2*cos(arg2)
      End Do

c---------------------------------
c specify the boundary values of f
c---------------------------------

      Do i=1,npts
       if(Ibc.eq.1) then
           g(i) = 0.0
       else if(Ibc.eq.2) then
           g(i) = p(i,3)
       else if(Ibc.eq.3) then
           g(i) = p(i,3)*p(i,3)
       end if
      End Do

c-----------------------------
c read the triangle quadrature
c-----------------------------

      call gauss_trgl
     +
     +  (mint
     +  ,xiq,etq,wq
     +  )

c----------------------------
c Compute alpha, beta, gamma,
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
c  xmon, ymom,zmom: x, y, and z surface moments
c                   over each element
c  area:            total surface area
c  vlm:             volume
c  crvmel:          mean curvature over each element
c  vna:             averaged normal vector at the nodes
c------------------------------------------------

      call elm_geom
     +
     +  (nelm,npts,mint
     +  ,xmom,ymom,zmom
     +  ,area,vlm
     +  ,cx,cy,cz
     +  )

c----
c normalize
c----

      arean = area/cell_area  

      write (6,*)
      write (6,110) arean
      write (6,111) vlmn
      write (6,*)

c-----------------
c printing session
c-----------------
c
c     write (6,*)
c     write (6,*)   " averaged normal vector:"
c     write (6,*)
c
c     Do i=1,npts
c      write (6,100) i,vna(i,1),vna(i,2),vna(i,3)
c     End Do
c-------------------------------------

c--------------------------
c prepare for the 3d_2p lgf
c--------------------------

      call ewald_3d_2p
     +
     +  (a11,a12,a21,a22
     +  ,b11,b12,b21,b22
     +  ,ew,area_lgf
     +  )

      write (6,*)
      write (6,*) "Ewald parameter: ",ew
      write (6,*) "Area of the unit cell: ",area_lgf
      write (6,*)

c---
c initialize the dipole
c to an arbitrary value
c---

      Do i=1,npts
       dpl(i) = 0.9876
c      dpl(i) = 0.1
      End Do

c--------------------------------------------
c iterative solution of the integral equation
c for the dipole strength
c--------------------------------------------

      Iter = 0

 71   Continue           ! target for iterations

      Iter = Iter+1

c-----------------------------------
c compute the double-layer potential
c at the nodes
c-----------------------------------

c     write (6,*) " Computing the principal value of the dlp"

      call ldlp_3d_2p
     +
     +  (npts
     +  ,nelm
     +  ,mint
     +  ,dpl
     +  ,dlp_pv
     +  )

c     write (6,*) " Done"

c---------------------
c update the dipole and compute the error
c---------------------

      error = 0.0D0

      Do i=1,npts

       save = dpl(i)

       dpl(i) = -cud*2.0*( dlp_pv(i) + grad_inf*p(i,3) - g(i) )

       error = error + (save-dpl(i))**2

c      write (6,100) i,g(i),dpl(i),dlp_pv(i)

      End Do

      error = sqrt(error)/npts

c----------------------------------
c compute the displacement constant
c----------------------------------

      Do i=1,npts
        auxs(i) = dpl(i)*vna(i,3)
      End Do

      call srf_int_3d_2p
     +
     +   (nelm
     +   ,npts
     +   ,mint
     +   ,auxs
     +   ,disp_c
     +   )

      disp_c = 0.5*disp_c/cell_area
c     disp_c = cud*0.5*disp_c/cell_area

      write (6,456) Iter,disp_c,error

c---
c more iterations?
c---

      if(error.le.error_max) Go to 711

      if(Iter.lt.Nter) Go to 71

c     write (6,*)
c     write (6,*) " One more iteration batch ? 1 for Yes, 0 for NO"
c     write (6,*) "-----------------------------------------------"
c     read  (5,*) more

      more = 0

      if(more.eq.1) then
       Iter = 0
       Go to 71
      end if

 711   Continue

c-------------------
c record the solution
c-------------------

      write (6,*)
      write (6,*)  " Solution converged"
      write (6,*)

c     write (6,*)
c     write (6,*) "converged solution: x-position, dipole"
c     write (6,*)

      Do i=1,npts
       write (4,100) i,p(i,1),dpl(i)
       write (6,100) i,p(i,1),dpl(i)
      End Do

c--
c record for visualization
c--

      open (1,file="ldr_3d_2p.net")

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

      Do k=1,nelm
       call printel (k,Index,dpl)  ! print in file "caps_3d.net"
      End Do

      close (1)

c----------------
c post-processing
c----------------

  97  Continue

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

c=======================
      if(menu.eq.1) then
c=======================

      write (6,*)
      write (6,*) " Enter the coordinates of the point"
      write (6,*) " where the solution is to be evaluated"
      write (6,*) "--------------------------------------"
      read  (5,*) xp,yp,zp

      call ldlpp_3d_2p
     +
     +  (npts
     +  ,nelm
     +  ,mint
     +  ,dpl
     +  ,xp,yp,zp
     +  ,ptl
     +  )

      write (6,*) "The solution is:"
      write (6,101) ptl

c===========
      end if
c===========

c=======================
      if(menu.eq.2) then
c=======================

      call sgrad_3d_2p
     +
     +   (npts,nelm
     +   ,f
     +   ,sgrad_f
     +   )

      write (6,*)
      write (6,*) "x-position - surface gradient"
      write (6,*)

      Do i=1,npts
       write (6,100) i,p(i,1),sgrad_f(i,1),sgrad_f(i,2),sgrad_f(i,3)
      End Do

c===========
      end if
c===========

c===========
      if(menu.eq.3) then
c===========

c---
c compute the strength of the
c vortex sheet
c from the dipole strength
c using the equation:
c
c zeta = n x grad(dpl)
c---

      call vs_3d_2p_circ_vort
     +
     +    (npts,nelm
     +    ,dpl
     +    ,zeta
     +    )

c---
c compute the vector potential
c at the nodes as a single-layer
c potential
c---

      call gauss_leg(NGL,zz,ww)

      call lslp3_3d_2p
     +
     +  (npts,nelm
     +  ,mint,NGL
     +  ,zeta
     +  ,slp3
     +  )

c---
c compute the normal component
c of the curl of slp3
c---

      call vs_3d_2p_curl
     +
     +   (npts,nelm
     +   ,slp3
     +   ,gradn_f
     +   )

c---
c add the basic field
c compute the projection of the curl of slp3
c onto the normal vector
c
c gradn_f_n = n.gradn_f
c---

      write (6,*)
      write (6,*) "x-position - magnitude of n.grad(f)"
      write (6,*)

      Do i=1,npts

       gradn_f(i,3) = grad_inf+gradn_f(i,3)   ! add the basic field

       gradn_f_n(i) = gradn_f(i,1)*vna(i,1)
     +               +gradn_f(i,2)*vna(i,2)
     +               +gradn_f(i,3)*vna(i,3)

c      write (6,100) i,p(i,1),gradn_f_n(i)

      End Do

c---
c compute the surface integral
c of gradn_f_n
c---

      call srf_int_3d_2p
     +
     +   (nelm,npts,mint
     +   ,gradn_f_n
     +   ,flow
     +   )

      write (6,*)
      write (6,457) flow
      write (6,*)

c---------------
c print elements
c--------------

      Do k=1,nelm

        write (3,105) nseven,k

        i = n(k,1)
        write (3,102) p(i,1),p(i,2),p(i,3),dpl(i),gradn_f_n(i)
        i = n(k,4)
        write (3,102) p(i,1),p(i,2),p(i,3),dpl(i),gradn_f_n(i)
        i = n(k,2)
        write (3,102) p(i,1),p(i,2),p(i,3),dpl(i),gradn_f_n(i)
        i = n(k,5)
        write (3,102) p(i,1),p(i,2),p(i,3),dpl(i),gradn_f_n(i)
        i = n(k,3)
        write (3,102) p(i,1),p(i,2),p(i,3),dpl(i),gradn_f_n(i)
        i = n(k,6)
        write (3,102) p(i,1),p(i,2),p(i,3),dpl(i),gradn_f_n(i)
        i = n(k,1)
        write (3,102) p(i,1),p(i,2),p(i,3),dpl(i),gradn_f_n(i)

      End Do

      write (3,100) null

c===========
      end if
c===========

      Go to 97    ! Return to the post-processing menu

c---
c done
c---

 99   Continue

      write (4,703) disp_c
      write (4,457) flow

      write (3,703) disp_c
      write (3,457) flow

      write (4,700) cud 
     +             ,grad_inf
     +             ,ndiv
     +             ,a11,a12
     +             ,a21,a22
     +             ,amp1,amp2
     +             ,mint
     +             ,NGL
     +             ,Method_lgf
     +             ,Max1,Max2,Max3
     +             ,Nter
     +             ,error_max

      write (3,700) cud
     +             ,grad_inf
     +             ,ndiv
     +             ,a11,a12
     +             ,a21,a22
     +             ,amp1,amp2
     +             ,mint
     +             ,NGL
     +             ,Method_lgf
     +             ,Max1,Max2,Max3
     +             ,Nter
     +             ,error_max

      close (2)
      close (3)

  100 Format (1x,i4,10(1x,f10.5))
  101 Format (f15.10)
  102 Format (10(1x,f10.5))
  105 Format (1x,i5,1x,i5)

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)

  139 Format (i4," : ",i1,";",9(1x,i4))

  456 Format ("Iteration: ",I3," c_d=",f12.8," Error: ",f20.15)
  457 Format ("Flow across the boundary:",f15.10)

  700 Format (/,
     +        "cud            =",f15.10,/,
     +        "grad_inf       =",f7.3,/,
     +        "Ndiv           =",i1,/,
     +        "a_1            =",f10.5,1x,f10.5,/,
     +        "a_2            =",f10.5,1x,f10.5,/,
     +        "amp1,amp2      =",f10.7,1x,f10.7,/,
     +        "mint           =",i2,/,
     +        "NGL            =",i2,/,
     +        "Method_lgf     =",i1,/,
     +        "Max1,Max2,Max3 =",i2,1x,i2,1x,i2,/,
     +        "Nter           =",i4,/,
     +        "Error          =",f15.12
     +        )

  703 Format ("Displacement constant",f15.10)

  800 Format(" centroid:",3(1x,f10.5))

      Stop
      End
