      program prtcl_3d_mob_dlr_se

c==========================================
c FEMLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------
c Mobility problem for a freely suspended particle
c computed by the double-layer represenation (dlr)
c
c The velocity of translation and angular velocity
c of rotation are computed by solving an integral
c equation of the second kind for the dipole strength
c distribution
c
c SYMBOLS:
c -------
c
c u(i,3):      incindent velocity at the ith node
c dpl(i,3):    dipole strength at the ith node
c dlp_pv(i,3): principal value of the double-layer potential
c              at the ith node
c ngsp:        number of global spectral nodes
c psp:         position of spectral nodes (global)
c xiq,etq,wq:  Gauss-triangle quadrature
c xir,etr,wr:  integration rule
c---------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  p(1026,3)
      Dimension ne(1026,7)
      Dimension  c(1026)
      Dimension  color(1026)

      Dimension n(512,6),nbe(512,3)

      Dimension  psp(20000,3)   ! spectral node
      Dimension  dpl(20000,3)

      Dimension vmaster(8)

      Dimension    van(500,500)
      Dimension vaninv(500,500)

      Dimension xiq(20),etq(20),wq(20)
      Dimension xir(250000),etr(250000),wr(250000)    ! integration weights

c---
c common blocks
c---

      common/points/p,ne
      common/elmnts/n,nbe

      common/var1/Iflow
      common/var2/Uinf,shrt
      common/var3/wall
      common/var4/Iflowinf

      common/spectral1/mpoly,npoly,vmaster,van,vaninv
      common/spectral2/nesp,ngsp
      common/spectral5/psp,dpl

      common/trq/xiq,etq,wq
      common/trr/xir,etr,wr
      common/trri/nbase

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
      spi4 = Dsqrt(pi4)

      null  = 0
      Nfour  = 4
      Nseven = 7
      oot  = 1.0D0/3.0D0

c-----------
c input data
c-----------

      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 to read data from file: prtcl_3d_mob_dlr_se.dat'
      write (6,*) ' 2 to type data into the keyboard'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*)   Ienrd

      if(Ienrd.eq.0) Go to 99

c------------------------
      if(Ienrd.eq.1) then
c------------------------

      open (2,file="prtcl_3d_mob_dlr_se.dat")

       read (2,*) Ioctaicos
       read (2,*) ndiv
       read (2,*)
       read (2,*) boa,coa
       read (2,*) req
       read (2,*) cxp
       read (2,*) cyp
       read (2,*) czp
       read (2,*) phi1,phi2,phi3
       read (2,*)
       read (2,*) Intm
       read (2,*) mint
       read (2,*)
       read (2,*) tol
       read (2,*) Nter
       read (2,*) 
       read (2,*) Iflow
       read (2,*) wall
       read (2,*) 
       read (2,*) Uinf
       read (2,*) shrt
       read (2,*) 
       read (2,*) mpoly
       read (2,*) Idist     ! 0 for uniform, 1 for spectral

      close (2)

c---------
      else
c---------

c-----------

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

      write (6,*) " Will discretize an ellipsoid"
      write (6,*) " into 6-node curved triangles"
      write (6,*)
      write (6,*) " Select the level of triangulation"
      write (6,*) " choose from 0, 1, 2, 3; 99 to quit"
      write (6,*) " ----------------------------------"
      read  (5,*) ndiv

      if(ndiv.eq.99) Go to 99

      if((ndiv.lt.0).or.(lt.gt.3))  then
        write (6,*) 'out of range; please try again'
        Go to 91
      end if

      write (6,*)
      write (6,*) " The x,y,z semi-axes of the ellipsoid are: a,b,c"
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
      write (6,*) "      for integrating over each triangle"
      write (6,*) " Enter m"
      write (6,*) " choose from 3, 4, 5, 6, 9, 12, 13"
      write (6,*) " Enter 0 to quit"
      write (6,*) "----------------"
      read  (5,*) mint

      write (6,*)
      write (6,*) " Enter the order of the quadrature"
      write (6,*) " for integration of the singular triangles"
      write (6,*) " -----------------------------------------"
      read  (5,*) NGL

      write (6,*)
      write (6,*) " Enter the tolerance for stopping"
      write (6,*) " the iterations"
      write (6,*) " --------------"
      read  (5,*) tol

      write (6,*)
      write (6,*) " How many batch iterations in solving"
      write (6,*) " the integral equation ?"
      write (6,*) "------------------------"
      read  (5,*) Nter

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 for flow in free space"
      write (6,*) " 2 for flow in a domain bounded by a plane wall"
      write (6,*) " ----------------------------------------------"
      read  (5,*) Iflow

      if(Iflow.eq.2) then

        write (6,*) " The wall is located at x=wall"
        write (6,*) " Please enter: wall"
        write (6,*) " ------------------"
        read  (5,*) wall

      end if

c-----------
      end if
c-----------

c-------------------------------------
c specify the incident velocity: u_inf
c-------------------------------------

   98 Continue

      write (6,*)
      write (6,*) 'Choose the ambient flow'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*)
      write (6,*) ' 1 for u_x = U'
      write (6,*) ' 2 for u_y = U'
      write (6,*) ' 3 for u_z = U'
      write (6,*) ' 4 for u_y = s x'
      write (6,*) ' 5 for u_z = s x'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) Iflowinf

      if(Iflowinf.eq.0) Go to 99

c------------------
c open output files
c------------------

      open (4,file="prtcl_3d_mob_dlr_se.out")

c-----------------------------------------
c read the triangle integration quadrature
c-----------------------------------------

      call gauss_trgl
     +
     +   (mint
     +   ,xiq,etq,wq
     +   )

c-----------------------------------
c read the triangle integration rule
c-----------------------------------

      if(Intm.eq.2) then
       call trgl_rule (mint,nbase,xir,etr,wr)
      end if

c----------------------------
c define the collocation nodes
c over the standard triangle in the xi-eta plane
c
c Compute the inverse
c of the vandermonde matrix
c---------------------------

       call vander (Idist)

c------------------
c solver the integral equation
c to compute the particle translational 
c and angular velocities
c------------------

      call prtcl_3d_mob_dlr_se_slv
     +
     +    (Ioctaicos
     +    ,ndiv
     +    ,boa,coa
     +    ,req
     +    ,cxp,cyp,czp
     +    ,phi1,phi2,phi3
     +    ,intm        ! integration method
     +    ,mint
     +    ,tol
     +    ,Nter
     + 
     +    ,npts,nelm   ! output ||
     +                 !        \/
     +    ,dpl_dn     
     +    ,Vx,Vy,Vz 
     +    ,Ox,Oy,Oz 
     +    )

c------
c print
c------

      write (6,*)

      if(Iflow.eq.1) then

       write (6,456) dpl_dn,Vx,Vy,Vz,Ox,Oy,Oz

      else if(Iflow.eq.2) then

       Oxr = 2.0D0*Ox   ! scale
       Oyr = 2.0D0*Oy
       Ozr = 2.0D0*Oz

       write (6,456) dpl_dn,Vx/cx,Vy/cx,Vz/cx,Oxr,Oyr,Ozr

      end if

c--------------------
c print all triangles
c--------------------

      open (1,file="prtcl_3d_mob_dlr_se.net")

      Index = 2
      Index = 1

      if(Index.eq.1) then          ! 6-node triangles
        write (1,*) Nseven
        write (1,*) 7*nelm
        write (1,*)   nelm
      else if(Index.eq.2) then     ! 3-node triangles
        write (1,*) Nfour
        write (1,*) 16*nelm
        write (1,*)  4*nelm
      end if

      Do i=1,npts
       color(i) = p(i,1)
      End Do

      Do k=1,nelm
       call printel (k,Index,color)  ! print in file 1
      End Do

      write (1,*) ngsp
c     write (6,*) ngsp

      Do i=1,ngsp
        write (1,102) psp(i,1),psp(i,2),psp(i,3)
c       write (6,102) psp(i,1),psp(i,2),psp(i,3)
      End Do

c     write (1,*) ntespn
c     Do k=1,nelm
c      Do i=1,nesp
c       write (1,102) xsp(k,i),ysp(k,i),zsp(k,i)
c      End Do
c     End Do
c     write (1,*) 2*nesp
c     Do i=1,nesp
c      write (1,102) xsp(1,i),ysp(1,i),zsp(1,i)
c     End Do
c     Do i=1,nesp
c      write (1,102) xsp(2,i),ysp(2,i),zsp(2,i)
c     End Do

      close (1)

c-------
c repeat
c-------

      Go to 98

c-----
c done
c-----

 99   Continue

      close (4)

  100 Format (1x,i5,10(1x,f10.5))
  101 Format (f15.10)
  102 Format (10(1x,f12.8))

  110 Format(" Surface Area :",F15.10)
  111 Format(" Volume       :",F15.10)

  456 Format (" Proj:",f12.8,/
     +       ," Vx: ",f12.8," Vy: ",f12.8," Vz: ",f12.8,/
     +       ," Ox: ",f12.8," Oy: ",f12.8," Oz: ",f12.8)

  800 Format(" centroid:",3(1x,f10.5))

      Stop
      End
