      program prtcl_2d

c===========================================
c FDLIB, BELMIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c-----------------------------------------------------------
c Two-dimensional Stokes flow
c past a collection of stationary, translating,
c or rotating  particles,
c for a variety of configurations.
c
c This program solves the integral equation for the traction
c over the particle contours,
c and then computes streamlines.
c
c
c SYMBOLS:
c -------
c
c Nprtcl: Number of particles
c
c NE(i):  Number of elements along the ith particle
c         contour
c
c Itp(i): Index for the elements on the ith particle:
c         1 for a straight segments
c         2 for native elements of circle or ellipses
c
c xcntr(i), ycnrt(i) : x and y coordinates
c                      of the center of the ith particle
c
c tilt(i) : Tilting angle of ith particle
c
c xw, yw:  coordinates of end-points of elements
c
c X0, Y0:  coordinates of collocation points
c t0:	   angle subtended from the center of a circular segment 
c          at collocation points
c S0:	   polygonal arc length at collocation points
c
c ux0,   uy0: perturbation velocity components at collocation points
c vnx0, vny0: normal vector at collocation points (into the flow)
c tnx0, tny0: tangent vector at collocation points (into the flow)
c
c Iprt(i):  Particle number hosting the ith collocation point
c
c elml(i):  Element arc length hosting the ith collocation point
c
c fx, fy(i,j): 	Cartesian components of the traction
c               on the jth element of the ith particle
c
c ftn(i,j):  Tangential component of the traction
c            on the jth element of the ith particle
c
c NGL: 	Number of Gaussian points for integration
c       over each element
c
c shrt:	  shear rate of incident flow
c delta:  pressure gradient of incident flow
c
c Nblocks:   Number of diagonal blocks in Master Linear System
c Lump(i):   Number of particle sub-blocks in the ith block
c
c forcex(i), forcey(i):	force on the ith particle
c torque(i):		torque on the ith particle
c
c Isolve:  1 solve according to prtcl_2d_sys
c            as described in the text
c
c          2 solve according to prtcl_2d_sys1
c            based on cloud particle clustering
c
c Iprec:  Singular preconditioning flag
c
c Ireg:  Regularization flag:
c
c        1: replace last equation for each particle
c           by the integral constranint Int f.n dl = 0
c        3: deflate the integral equation
c
c Ireduce: superseded by Ireg
c
c
c Capacity:
c ---------
c
c       72 particles
c       64 elements along each particle
c
c---------------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(72),  Itp(72)
      Dimension axis1(72),axis2(72)
      Dimension xcntr(72),ycntr(72),tilt(72)

      Dimension  xw(72,65),yw(72,65),tw(72,65)
      Dimension  sw(72,65)

      Dimension   X0(4608),  Y0(4608),T0(4608)
      Dimension   S0(4608)
      Dimension  ux0(4608), uy0(4608)
      Dimension elml(4608)
      Dimension vnX0(4608),vnY0(4608)
      Dimension tnX0(4608),tnY0(4608)
      Dimension Iprt(4608)

      Dimension  fx(72,64),fy(72,64),ftn(72,64)

      Dimension forcex(72),forcey(72),torque(72)

c---------------------------------- 
c for the master linear system (MLS)
c---------------------------------- 

      Dimension Lump(72)

c--------------------------
c Gauss-Legendre quadrature
c--------------------------

      Dimension ZZ(20),WW(20)

c---
c streamlines
c---

      Dimension xstr(900),ystr(900)   ! for streamlines

c--- 
c for sgf_2d_2p
c--- 

      Dimension  qxx(-15:15,-15:15)
      Dimension  qxy(-15:15,-15:15)
      Dimension  qyy(-15:15,-15:15)
      Dimension vxxx(-15:15,-15:15)
      Dimension vxxy(-15:15,-15:15)
      Dimension vyxy(-15:15,-15:15)
      Dimension vxyx(-15:15,-15:15)
      Dimension vxyy(-15:15,-15:15)
      Dimension vyyy(-15:15,-15:15)
      Dimension  ppx(-15:15,-15:15)
      Dimension  ppy(-15:15,-15:15)

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NE,Itp,Ncl
      common/INTGR2/NGL

      common/REAL1/visc,shrt,delta,Uinf,wall
      common/REAL2/Uprtcl,Vprtcl,Aprtcl

      common/CHANNELI/IQPD,NGFWW
      common/CHANNELR/U1,U2,RL,h

      common/pax/axis1,axis2
      common/pap/xcntr,ycntr,tilt

      common/points/xw,yw,tw

      common/colloc1/ux0,uy0
      common/colloc2/vnx0,vny0,elml
      common/colloc3/tnx0,tny0
      common/colloc4/x0,y0,t0
      common/colloc5/Iprt
      common/colloc6/fx,fy

c---
c for sgf_2d_2p
c---
      common/aaaa/a11,a12,a21,a22
      common/bbbb/b11,b12,b21,b22
      common/ewew/ew,tau
      common/mmmm/max1,max2
      common/qqqq_2d/qxx,qxy,qyy
      common/vvvv_2d/vxxx,vxxy,vyxy,vxyx,vxyy,vyyy,ppx,ppy
c---
c various
c---
      common/ZZWW/ZZ,WW
      common/piii/pi,pi2,pi4

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi
      pi6 = 6.0D0*pi
      pi8 = 8.0D0*pi

      Null = 0
      None = 1
      Ntwo = 2

      zero  = 0.0D0
      five  = 5.0D0
      fivem =-5.0D0

c------------
c preferences
c------------

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 to solve the integral equation"
      write (6,*) " 2 to draw streamlines"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Menu

      If(Menu.eq.0) Go to 99

  94  Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) 
      write (6,*) " 1 for shear flow past a swarm of particles "
      write (6,*) "   translating above a wall"
      write (6,*) 
      write (6,*) " 2 for shear flow past a periodic array "
      write (6,*) "   of particles translating above a wall"
      write (6,*) 
      write (6,*) " 3 for flow past a periodic array "
      write (6,*) "   of particles translating in a channel"
      write (6,*) 
      write (6,*) "10 for flow past a doubly-periodic array"
      write (6,*) "----------------------------------------"
      write (6,*) 
      read  (5,*) Iflow

      If(Iflow.eq.0) Go to 99

      If(Iflow.ne.1.and.Iflow.ne.2
     +             .and.Iflow.ne.3
     +             .and.Iflow.ne.10
     +  ) then
       write (6,*)
       write (6,*) " Invalid selection; please repeat"
       write (6,*)
       Go to 94
      End If

c-----------------
c  Read parameters
c-----------------

      open (4,file='prtcl_2d.dat')

        read (4,*) visc     ! fluid viscosity
        read (4,*) 
        read (4,*) NGL      ! Gauss Legendre integration points
        read (4,*) 
        read (4,*) shrt     ! shear rate
        read (4,*) delta    ! pressure gradient
        read (4,*) wall     ! the wall is located at y = wall
        read (4,*) RL       ! period
        read (4,*) 
        read (4,*) Uprtcl   ! Particle translational x-velocity
        read (4,*) Vprtcl   ! Particle translational y-velocity
        read (4,*) Aprtcl   ! Particle angular velocity about its axis
        read (4,*) 
        read (4,*) U1       ! channel flow: lower wall velocity
        read (4,*) U2       ! channel flow: upper wall velocity
        read (4,*) h        ! channel flow: channel semi-width
        read (4,*) IQPD     ! specify constant pressure drop or flow rte
        read (4,*) NGFWW    ! runcation limit for sgf_2d_ww
        read (4,*) 
        read (4,*) Uinf       ! velocity if incident flow
        read (4,*) 
        read (4,*) a11,a12    ! first base vector
        read (4,*) a21,a22    ! second base vector
        read (4,*) max1,max2  ! truncation limits for the Green's function
        read (4,*) 
        read (4,*) Isolve ! solution method (1 or 2)
        read (4,*) 
        read (4,*) Method ! choose the integral formulation (for Isolve=1)
        read (4,*) 
        read (4,*) Iprec  ! preconditioning flag (for Isolve=1)
        read (4,*) Ireg   ! regularization flag  (for Isolve=1) (replaced Ireduce)
        read (4,*) 
        read (4,*) expn  ! expansion factor for clusters (for Isolve=2)
        read (4,*) 
        read (4,*) Niter  ! max number of iterations
        read (4,*) tol    ! accuracy

      close (4)

c--------------------------------------
c  Read particle geometrical properties
c--------------------------------------

      If(Iflow.eq.1)  open (4,file='prtcl_2d_geo01.dat')
      If(Iflow.eq.2)  open (4,file='prtcl_2d_geo02.dat')
      If(Iflow.eq.3)  open (4,file='prtcl_2d_geo03.dat')
      If(Iflow.eq.10) open (4,file='prtcl_2d_geo10.dat')

      read (4,*) Nprtcl       ! number of particles

      Do i=1,Nprtcl
         read (4,*) idle,axis1(i),axis2(i)
     +                  ,xcntr(i),ycntr(i),tilt(i)
     +                  ,NE(i),Itp(i)
      End Do

c------------
c comment out one of the following two options
c-----------

c--- option 1: hard reset

      Nblocks = Nprtcl
      Do i=1,Nblocks
        Lump(i) = 1
      End Do

c--- option 2: read

c     read (4,*) Nblocks               ! number of particle blocks
c     read (4,*) (Lump(i),i=1,Nblocks) ! number of particles within a block

      close (4)

c-----------------
c open output file
c-----------------

      open (3,file="prtcl_2d.xy")

c--------
c prepare
c--------

      Do i=1,Nprtcl
       tilt(i) = tilt(i)*pi     ! particle tilting angles
      End Do

c-----------------
c call quadratures
c-----------------

      call gauss_legendre (NGL,ZZ,WW)

c-------------
c channel flow
c-------------

      if(Iflow.eq.3)  then
        h2 = 2.0D0*h 
      End If

c-----------------------------------
c plot walls in file: prtcl_2d.strml
c-----------------------------------

      open (7,file="prtcl_2d.strml")

      If(Iflow.eq.1.or.Iflow.eq.2) then  ! plane wall

        write (7,*)   Ntwo
        write (7,104) None,fivem,wall
        write (7,104) Ntwo,five ,wall

      Else If(Iflow.eq.3) then           ! two channel walls

        write (7,*)   Ntwo
        write (7,104) None,fivem,-h2
        write (7,104) Ntwo,five , h2
        write (7,*)   Ntwo
        write (7,104) None,fivem,-h2
        write (7,104) Ntwo,five , h2

      End If

c---------------------------------
c Generate element end-points
c and record in file: prtcl_2d.str
c---------------------------------

      write (6,*) " prtcl_2d: Generating element end-points"

      Do i=1,Nprtcl

        cs = Dcos(tilt(i))
        sn = Dsin(tilt(i))

        Dtt = pi2/NE(i)

        Do j=1,NE(i)+1
         tt      = (j-1.0D0)*Dtt
         TW(i,j) = tt
         tmpx    = axis1(i)*Dcos(tt)
         tmpy    = axis2(i)*Dsin(tt)
         xw(i,j) = xcntr(i)+tmpx*cs-tmpy*sn     ! rotate to tilt
         yw(i,j) = ycntr(i)+tmpx*sn+tmpy*cs
        End Do

        write (7,*) NE(i)+1,i

        Do j=1,NE(i)+1
         write (7,104) j,xw(i,j),yw(i,j)
        End Do

      End Do

c------------------------------
c Count the collocation points;
c one per element
c------------------------------

      Ncl = 0

      Do i=1,Nprtcl
       Ncl = Ncl+NE(i)
      End Do

      write (6,*) " prtcl_2d:" ,Ncl," collocation points"

c---------------------------------------------------
c  Compute the polygonal arc length at end-points
c  and element collocation points (mid-points)
c  along each particle contour;
c  will be used for plotting purposes
c---------------------------------------------------

      Ic = 0     ! collocation point counter

      Do i=1,Nprtcl

        sw(i,1) = 0.0D0   ! arc length at first end-point

        Do j=1,NE(i)

         Ic = Ic+1
         sw(i,j+1) = sw(i,j)
     +             +  Dsqrt( (xw(i,j+1)-xw(i,j))**2
     +                      +(yw(i,j+1)-yw(i,j))**2 )
         S0(Ic) = 0.5D0*(sw(i,j+1)+sw(i,j))

        End Do

      End Do

c----------------------------------------
c Define and confirm solution block sizes
c----------------------------------------

      If(Nblocks.eq.1) then           ! only one block

        Lump(1) = Nprtcl              ! lump all particles

      Else If(Nblocks.eq.Nprtcl) then ! each particle is a block

        Do i=1,Nprtcl                 ! each particle is an individual
         Lump(i) = 1
        End Do

      Else

        Nsum = 0                        ! check for consistency
        Do i=1,Nblocks
         Nsum = Nsum + Lump(i)
        End Do

        If(Nsum.ne.Nprtcl) then 
          write (6,*)
          write (6,*) " prtcl_2d: Inappropriate partitioning"
          write (6,*) "           of the Master Linear System."
          write (6,*) "           Quitting."
          write (6,*)
          stop
        End If

      End If

c-------------------------------------------
c Collocation points:
c
c Compute: coordinates, normal and tangential vector,
c          element lengths
c          host particle of collocation points
c          disturbance velociy
c-------------------------------------------

      call coll_points (Iflow)

c-------------------------------
c Display the collocation points
c-------------------------------
c
c     Ic = 0
c     Do i=1,Nprtcl
c       write (3,*) NE(i)
c       write (6,*) NE(i)
c       Do j=1,NE(i)
c         Ic = Ic+1
c         write (6,102) j,X0(Ic),Y0(Ic),S0(Ic)
c         write (3,102) j,X0(Ic),Y0(Ic),S0(Ic)
c        End Do
c       End Do
c     End Do
c-----------------------------

c---------------------------------
c prepare for doubly periodic flow
c---------------------------------

      if(Iflow.eq.10) then

       call sgf_2d_2p_ewald
     +
     +     (a11,a12,a21,a22
     +     ,b11,b12,b21,b22
     +     ,ew,tau
     +     )

       call sgf_2d_2p_qqq
     +
     +     (b11,b12
     +     ,b21,b22
     +     ,max2,ew
     +     )

c---
       if(method.eq.2) then

        call sgf_2d_2p_vvv
     +
     +    (b11,b12
     +    ,b21,b22
     +    ,max2,ew
     +    )
c---

       End If

      End If

c--------------------------------
c Prepare to draw streamlines
c
c Read collocation point position
c and tractions from input file:
c
c prtcl_2d.inp
c--------------------------------

      if(menu.eq.2) then

        open (8,file="prtcl_2d.inp")

        Ic = 0

        Do i=1,Nprtcl

          read (8,*) NE(i)

          Do J=1,NE(i)
           Ic = Ic+1
           read (8,102) Idle,X0(Ic),Y0(Ic),S0(Ic)
     +                      ,fx(I,J),fy(I,J),ftn(I,J)
          End Do
        End Do

        close (8)

        Go to 92    ! proceed to compute streamlines

      End If

c---------
c  Warning
c---------

      if(Iprec.eq.1.and.Ireg.eq.0) then
        write (6,*)
        write (6,*) " prtcl_2d: WARNING"
        write (6,*)
        write (6,*) " Singular preconditioning without regularization"
        write (6,*) " selected; it very likely that the computation"
        write (6,*) " will fail"
c       pause
      end If

c--------------------------------
c  Generate and solve
c  the Master Linear System (MLS)
c--------------------------------

      if(Isolve.eq.1) then

      call prtcl_2d_sys
     +
     +    (Iflow
     +    ,Method
     +    ,Iprec
     +    ,Ireg
     +    ,Nblocks
     +    ,Lump
     +    ,Niter
     +    ,tol
     +    ,Istop
     +    )    

      else If(Isolve.eq.2) then

       call prtcl_2d_sys1
     +
     +   (Iflow
     +   ,expn
     +   ,Niter
     +   ,tol
     +   ,Istop
     +   )

      else

         write (6,*) " prtcl_2d: solution method not selected"
         stop

      end if

c--------------------------------------
c Compute the total traction by adding
c         the traction of the incident flow
c
c This is done only for Method = 2,
c for then we solve the integral equation
c for the ``disturbance'' traction
c----------------------------------

      if(method.eq.2) then 

      Ic = 0

      Do I=1,Nprtcl
       Do J=1,NE(I)

       Ic = Ic+1

       If(Iflow.eq.1.or.Iflow.eq.2
     +   ) then

         yref = Y0(Ic)-wall
         sg_xx = -delta*X0(Ic)         ! stress tensor 
         sg_xy = visc*shrt+delta*yref  ! of the reference flow
         sg_yx = sg_xy
         sg_yy = sg_xx

       else if(Iflow.eq.3) then

         sg_xx = -delta*X0(Ic)                 ! stress tensor 
         sg_xy = visc*(U2-U1)/h2+delta*Y0(Ic)  ! of the reference flow
         sg_yx = sg_xy
         sg_yy = sg_xx

       else if(Iflow.eq.10) then

         sg_xx = 0.0D0
         sg_xy = 0.0D0
         sg_yx = 0.0D0
         sg_yy = 0.0D0

       end if

       fx(i,j) = fx(i,j)+sg_xx*vnx0(Ic)+sg_xy*vny0(Ic)
       fy(i,j) = fy(i,j)+sg_xy*vnx0(Ic)+sg_xx*vny0(Ic)

       end do
      end do

      end if

c----------------------------------------------
c compute:
c
c  (a) the force and torque on each particle
c  (b) the tangential components of the traction
c----------------------------------------------

      Ic = 0   ! collocation point counter

      Do i=1,Nprtcl

       forcex(i) = 0.0D0
       forcey(i) = 0.0D0
       torque(i) = 0.0D0

       Do j=1,NE(I)

        Ic = Ic+1

        forcex(i) = forcex(i) + fx(i,j)*elml(Ic)
        forcey(i) = forcey(i) + fy(i,j)*elml(Ic)

        torque(i) = torque(i)
     +            + ( (X0(Ic)-xcntr(i))*fy(i,j)
     +               -(Y0(Ic)-ycntr(i))*fx(i,j)
     +              )*elml(Ic)

        ftn(i,j) = fx(i,j)*tnx0(Ic)+fy(i,j)*tny0(Ic)

       End Do

      End Do

c-----------------
c Printing session
c-----------------

      write (6,*)  
      write (6,*)  " prtcl_2d: x, y, l, tractions, shear stress"
      write (6,*) 

      Ic = 0           ! collocation point counter

      Do i=1,Nprtcl

        write (6,*) NE(i),i
        write (3,*) NE(i),i

        Do J=1,NE(i)

         Ic = Ic+1

         sprint = S0(Ic)/sw(i,NE(i)+1)

         write (6,102) J,X0(Ic),Y0(Ic),sprint
     +                  ,fx(I,J),fy(I,J)
     +                  ,ftn(I,J)

         write (3,102) J,X0(Ic),Y0(Ic),sprint
     +                  ,fx(I,J),fy(I,J)
     +                  ,ftn(I,J)

        End Do
      End Do

      write (3,*) Null
      write (3,*)
      write (3,*)  " Particle No, Force and Torque"
      write (3,*)
      write (6,*)
      write (6,*)  " Particle No, Force and Torque"
      write (6,*)

      write (3,*) Nprtcl

      Do i=1,Nprtcl
       write (6,102) i,forcex(i),forcey(i),torque(i)
       write (3,102) i,forcex(i),forcey(i),torque(i)
      End Do

      write (3,*) Null
      write (3,*) 

c-------------------------------------
c END OF SOLVING THE INTEGRAL EQUATION
c-------------------------------------

      write (6,*)
      write (6,*) " Enter 1 to draw streamlines"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      Read  (5,*) menu

      If(menu.eq.0) Go to 99

  92  Continue

      write (6,*)
      write (6,*) " Enter the maximum number of points along "
      write (6,*) " a streamline before inquiring "
      write (6,*) " ------------------------------"
      read  (5,*) Mstr

      write (6,*)
      write (6,*) " Enter: "
      write (6,*)
      write (6,*) " 2 to select the second-order RK method"
      write (6,*) " 4 to select the fourth-order RK method"
      write (6,*) " 0 to quit"
      write (6,*) " -------- "
      read  (5,*) IRK

      If(IRK.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter the spatial step"
      write (6,*) " ----------------------"
      read  (5,*) spat_step

      open (9,file="strml.dat")

  22  Continue

      read (9,*) Xstart,Ystart    ! will stop if Xstart = 99 
                                  !           or Ystart = 99

      If(abs(Xstart-99).lt.0.0000001) Go to 99
      If(abs(Ystart-99).lt.0.0000001) Go to 99

      call prtcl_2d_strml
     +
     +     (Iflow
     +     ,Xstart,Ystart
     +     ,Mstr
     +     ,IRK
     +     ,spat_step
     +     ,Nsp
     +     ,xstr,ystr
     +     )

      If(Nsp.eq.0) Go to 22      ! return for another streamline

c-------------------------
c Print out the streamline
c-------------------------

      write (6,*)
      write (6,*) " One streamline with ",Nsp," points completed"
      write (6,*)

      write (7,104) Nsp   ! number of points along the streamline

      Do i=1,Nsp
        write (7,104) i,xstr(i),ystr(i)
      End Do

c-------------
      Go to 22       ! return for another streamline
c-------------

  99  Continue

c-----
c Done
c-----

      write (3,*) Null
      write (7,*) Null

c---
c record parameters
c---

      write (3,169) Iflow
      write (3,170) visc
      write (3,171) NGL
      write (3,172) shrt
      write (3,173) delta
      write (3,174) wall
      write (3,175) RL
      write (3,176) U1
      write (3,177) U2
      write (3,178) h
      write (3,179) IQPD
      write (3,180) NGFWW
      write (3,181) Uinf
      write (3,182) a11,a12
      write (3,183) a21,a22
      write (3,184) max1,max2
      write (3,185) Method
      write (3,186) Iprec
      write (3,986) Ireg
      write (3,188) Nblocks
      write (3,189) (Lump(i),i=1,Nblocks)
      write (3,190) Niter
      write (3,990) tol
      write (3,192) Nprtcl

      write (3,*) 
      write (3,*)  " i, xc, yc, axis1, axis2, tilt, NE, Itp"

      Do i=1,Nprtcl
       write (3,193) i,axis1(i),axis2(i)
     +                ,xcntr(i),ycntr(i),tilt(i)
     +                ,NE(i),Itp(i)
      End Do

c------------
c close files
c------------

      close (3)
      close (7)
      close (9)

c-----
c Done
c-----

 100  Format (1x,i3,1x,i3,10(f10.5))

 101  Format (80(1x,f7.3))
 102  Format (1x,i3,20(1x,f10.5))
 103  Format (1x,i3,5(1x,f15.10))
 104  Format (1x,i3,20(1x,f9.5))
 111  Format (80(10(1x,f7.3),/))
 112  Format (80(1x,f10.7))
 150  Format (1X,10(1X,f10.5))

 169  Format (1x,"Iflow     = ",I3)
 170  Format (1x,"Viscosity = ",f10.5)
 171  Format (1x,"NGL       = ",I3)
 172  Format (1x,"Shear rate= ",f10.5)
 173  Format (1x,"Delta     = ",f10.5)
 174  Format (1x,"Wall      = ",f10.5)
 175  Format (1x,"RL        = ",f10.5)
 176  Format (1x,"U1        = ",f10.5)
 177  Format (1x,"U2        = ",f10.5)
 178  Format (1x,"h         = ",f10.5)
 179  Format (1x,"IQPD      = ",I3)
 180  Format (1x,"NGFWW     = ",I3)
 181  Format (1x,"Uinf      = ",f10.5)
 182  Format (1x,"a11,a12   = ",f10.5,1x,f10.5)
 183  Format (1x,"a21,a22   = ",f10.5,1x,f10.5)
 184  Format (1x,"max1, max2= ",i2,1x,i2)
 185  Format (1x,"Method    = ",i1)
 186  Format (1x,"Iprec     = ",i1)
 986  Format (1x,"Ireg      = ",i1)
 188  Format (1x,"Nblocks   = ",i3)
 189  Format (1x,"lump      = ",50(1x,i2))
 190  Format (1x,"Niter     = ",i5)
 990  Format (1x,"Tol       = ",f15.10)
 192  Format (1x,"Nprtcl    = ",i5)
 193  Format (1x,i2,1x,5(1x,f10.5),1x,2(1x,i2))

 200  Format (1x,i3,1x,i3,10(f10.5))
 201  Format (1x,i3,10(f10.5))

 888  Format (1x,"Iter: ",i3," Point corr: ",f15.10)

      Stop
      End
