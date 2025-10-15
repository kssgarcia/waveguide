      program prtcl_ax

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c-----------------------------------------------------------
c Axisymmetric Stokes flow
c past a collection of stationary or translating particles
c for a variety of configuration
c
c This program solves the integral equation for the traction
c and then computes streamlines
c
c
c SYMBOLS:
c --------
c
c Nprtcl: Number of particles
c NE(i):  Number of elements along the trace of the ith particle
c
c Itp(i): Index for the elements on the ith particle:
c         1 for a straight segments
c         2 for native elements of ellipses 
c
c Isi(i): Shape index for the ith particle:
c         0 for compact particles with singly-connected contours
c         1 for particles with multiply-connected contours
c
c (xw, yw)   coordinates of end-points of elements
c
c (X0, Y0):  coordinates of collocation points
c t0:	     angle subtended from the center of a circular segment 
c S0:	     polygonal arc length
c
c (ux0, uy0):  perturbation velocity components at collocation points
c (vnx0, vny0): normal vector at collocation points (into the flow)
c (tnx0, tny0): tangent vector at collocation points (into the flow)
c
c arel(i): Area of axisymmetric surface hosting i collocation point
c
c fx, fy(,j):  Cartesian components of the traction
c              on the jth element of the ith particle.
c
c ftn: 	Tangential component of the traction
c       on the jth element of the ith particle.
c
c NGL: 	Number of Gaussian points for integration over each element
c
c Nblocks: Number of diagonal blocks in Master Linear System (MLS)
c Lump(i): Number of particle sub-blocks in the ith block
c Idmn(i):  Low  index of ith particle in master matrix
c Idmx(i):  High index of ith particle in master matrix
c
c forcex(i): force on the ith particle
c
c Capacity:
c --------
c
c    25 particles
c    64 elements along each particle contour
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(25),  Itp(25),Isi(25)
      Dimension xcntr(25),ycntr(25)
      Dimension axis1(25),axis2(25),tilt(25)

      Dimension xw(25,65),yw(25,65),tw(25,65)
      Dimension  sw(25,65)

      Dimension   X0(3200),  Y0(3200),T0(3200)
      Dimension S0(0:3200)
      Dimension  ux0(3200), uy0(3200)

      Dimension elar(3200)
      Dimension vnX0(3200),vnY0(3200)
      Dimension tnX0(3200),tnY0(3200)
      Dimension Iprt(3200)

      Dimension fx(25,64),fy(25,64)
      Dimension ftn(25,64)

      Dimension forcex(25)

c--- 
c for the master linear system (MLS)
c--- 

      Dimension Lump(25)

c---
c various
c---

      Dimension ZZ(20),WW(20)

c---
c streamlines
c---

      Dimension xstr(900),ystr(900)    ! for streamlines

c--------------
c common blocks
c--------------

      common/INTGR1/Nprtcl,NGL,NE,Itp,Ncl

      common/points/xw,yw,tw
      common/particles/xcntr,ycntr,axis1,axis2,tilt

      common/colloc1/fx,fy,ux0,uy0
      common/colloc2/vnx0,vny0,elar
      common/colloc3/tnx0,tny0
      common/colloc4/x0,y0,t0
      common/colloc5/Iprt

      common/REAL1/visc,Uinf,wall,pg,RL,sc
      common/REAL2/Uprtcl

      common/stream/xstr,ystr        ! for streamlines

c---
c for sgfaxct
c---

      common/sgfaxct/Nsum,Np

c---
c quadratures
c---

      common/ZZWW/ZZ,WW

c---
c various
c---

      common/piii/pi,pi2,pi4,pi8

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
      write (6,*) " Enter: "
      write (6,*) 
      write (6,*) " 1 to solve the integral equation"
      write (6,*) " 2 to draw streamlines"
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

  83  write (6,*) 
      write (6,*) " SELECT THE FLOW"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for uniform flow past stationary"
      write (6,*) "   or translating particles" 
      write (6,*) 
      write (6,*) " 3 for flow due to axial translation"
      write (6,*) "   normal to a plane wall"
      write (6,*) 
      write (6,*) " 5 for uniform flow past a translating file"
      write (6,*) "   of particles along the axis of a circular tube"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) Iflow

      if(Iflow.eq.0) Go to 99

c---
c trap
c---

      if(Iflow.ne.1.and.Iflow.ne.3
     +             .and.Iflow.ne.5
     +   ) then
       write (6,*) " prtcl_ax: Unfortunate selection"
       write (6,*) " prtcl_ax: Please try again"
       Go to 83
      end if

c--------------------------------------------------

c----------------
c read parameters
c----------------

      open (4,file='prtcl_ax.dat')

        read (4,*) visc      ! fluid viscosity
        read (4,*) 
        read (4,*) NGL       ! Gauss Legendre integration points
        read (4,*) 
        read (4,*) Uprtcl    ! particle translational velocity
        read (4,*) Uinf      ! velocity of incident flow
        read (4,*) 
        read (4,*) wall      ! the wall is located at x = wall
        read (4,*) 
        read (4,*) pg        ! pressure drop
        read (4,*) sc        ! tube radius
        read (4,*) RL        ! period
        read (4,*) Nsum      ! for computing sgf_ax_ct
        read (4,*) Np        ! for computing sgf_ax_ct
        read (4,*) 
        read (4,*) Iprec     ! preconditioning flag
        read (4,*) Ireduce   ! reduction flag
        read (4,*) 
        read (4,*) Niter                   ! max number of iterations
        read (4,*) tol                     ! accuracy

      close (4)

c-------------------------
c read particle properties
c-------------------------

      open (4,file='particle.dat')

        read (4,*) Nprtcl     ! number of particles

        Do i=1,Nprtcl
         read (4,*) idle,xcntr(i),ycntr(i)
     +                  ,axis1(i),axis2(i),tilt(i)
     +                  ,NE(i),Itp(i),Isi(i)
        End Do

        read (4,*) Nblocks               ! number particle blocks
        read (4,*) (Lump(i),i=1,Nblocks) ! number of particles inside a block

      close (4)

c------------------
c open output files
c------------------

      open (3,file="prtcl_ax.xy")
      open (7,file="prtcl_ax.str")

c--------
c prepare
c
c cross-section tilting angles
c for toroidal particles
c--------

      Do i=1,Nprtcl
        tilt(i) = tilt(i)*pi     
      End Do

c---
c quadrature
c---

      call gauss_leg (NGL,ZZ,WW)

c--------------------------------------------
c plot walls in streamline file: prtcl_ax.str
c--------------------------------------------

      if(Iflow.eq.3) then              ! trace of the wall

       write (7,*)   Ntwo
       write (7,104) None,wall,zero
       write (7,104) Ntwo,wall,five

      elseif(Iflow.eq.5) then         ! trace of the tube

       write (7,*)   Ntwo
       write (7,104) None,fivem,sc
       write (7,104) Ntwo,five ,sc

      end if
    
c-------------------------------------
c Generate particle-element end points 
c-------------------------------------

      Do i=1,Nprtcl    ! run over particles

        if(Isi(i).eq.0) then     ! singly connected shape (spheroid)
           tilt(i)  = 0.0D0
           ycntr(i) = 0.0D0
           Dtt = pi/NE(i)
           cs  = 1.0D0
           sn  = 0.0D0
        else                     ! multiply connected shape (torus)
           Dtt = pi2/NE(i)
           cs  = Dcos(tilt(i))
           sn  = Dsin(tilt(i))
        end if

        Do j=1,NE(i)+1
         Ic      = Ic+1
         tt      = (j-1.0D0)*Dtt
         TW(i,j) = tt
         tmpx    = axis1(i)*Dcos(tt)
         tmpy    = axis2(i)*Dsin(tt)
         XW(i,j) = xcntr(i)+tmpx*cs-tmpy*sn     ! rotate to tilt
         YW(i,j) = ycntr(i)+tmpx*sn+tmpy*cs
        End Do

      End Do

c----------------------------------------
c record in streamline file: prtcl_ax.str
c----------------------------------------

      Do i=1,Nprtcl

        write (7,*) NE(i)+1,i

        Do j=1,NE(i)+1
         write (7,104) j,xw(i,j),yw(i,j)
        End Do

        write (7,*) NE(i)+1,i    ! image contour
        Do j=1,NE(i)+1
         write (7,104) j,xw(i,j),-yw(i,j)
        End Do

      End Do

c------------------------------
c Count the collocation points
c one per element
c------------------------------

      Ncl = 0

      Do i=1,Nprtcl
       Ncl = Ncl+NE(i)
      End Do

      write (6,*) "prtcl_ax: ",Ncl," collocation points"

c--------------------------------------------
c Compute the polygonal arc length at element
c end  and mid points around each particle contour
c to be used for plotting purposes
c--------------------------------------------

      Ic = 0     ! collocation point counter

      Do i=1,Nprtcl

        sw(i,1) = 0.0D0    ! arc length at first end-point

        Do j=1,NE(i)

         Ic = Ic+1
         sw(i,j+1) = sw(i,j)+  Dsqrt((xw(i,j+1)-xw(i,j))**2
     +                              +(yw(i,j+1)-yw(i,j))**2)

         S0(Ic) = 0.5D0*(sw(i,j+1)+sw(i,j))   ! mid point

        End Do

      End Do

c----------------------------------------
c Define and confirm solution block sizes
c----------------------------------------

      if(Nblocks.eq.1) then             ! only one block

        Lump(1) = Nprtcl                ! lump all particles

      elseif(Nblocks.eq.Nprtcl) then   ! each particle is a block

        Do i=1,Nprtcl                   ! each particle is an individual
         Lump(i) = 1
        End Do

      else

        Nsum = 0                    ! check for consistency
        Do i=1,Nblocks
         Nsum = Nsum + Lump(i)
        End Do

        if(Nsum.ne.Nprtcl) then
          write (6,*) " prtcl_ax: inappropriate partitioning"
          write (6,*) "           of the Master Linear System"
          stop
        end if

      end if

c-------------------------------------------
c Collocation points:
c
c Compute:
c
c coordinates, normal and tangential vector,
c element lengths, and area of the surface of
c revolution
c-------------------------------------------

      call prtcl_ax_geo (Iflow)

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

c--------------------------------
c Prepare to draw streamlines
c
c Read collocation point position
c and tractions
c--------------------------------

      if(menu.eq.2) then

        open (8,file="prtcl_ax.inp")

        Do i=1,Nprtcl

         read (8,*) NE(i)

          Do J=1,NE(i)
           Ic = Ic+1
           read (8,102) Idle,X0(Ic),Y0(Ic),S0(Ic)
     +                      ,fx(I,J),fy(I,J)
     +                      ,ftn(I,J)
          End Do
        End Do

        close (8)

        Go to 92       ! proceed to compute streamlines

      end if

c--------------------------------
c  Generate and solve
c  the Master Linear System (MLS)
c--------------------------------

      call prtcl_ax_sys
     +
     +    (Iflow
     +    ,Iprec
     +    ,Ireduce
     +    ,Nblocks
     +    ,Lump
     +    ,Niter
     +    ,tol
     +    ,Istop
     +    )

c--------------------------------------
c tangential components of the traction
c and x component of the force exerted
c on each particle
c--------------------------------------

      Ic = 0

      Do I=1,Nprtcl

       forcex(i) = 0.0D0

       Do J=1,NE(I)
        Ic = Ic+1
        ftn(i,j)  = fx(i,j)*tnx0(Ic)+fy(i,j)*tny0(Ic)
        forcex(i) = forcex(i) + fx(i,j)*elar(Ic)
       End Do

      End Do

c-----------------
c Printing session
c-----------------

      write (6,*)  
      write (6,*)  " prtcl_ax: x, y, s, tractions, shear stress:"
      write (6,*) 

      Ic = 0           ! counter

      Do i=1,Nprtcl

        write (6,*) NE(i),i
        write (3,*) NE(i),i

        Do J=1,NE(i)
         Ic = Ic+1
         write (6,102) J,X0(Ic),Y0(Ic),S0(Ic)
     +         ,fx(I,J),fy(I,J)
     +         ,ftn(I,J)
         write (3,102) J,X0(Ic),Y0(Ic),S0(Ic)
     +         ,fx(I,J),fy(I,J)
     +         ,ftn(I,J)
        End Do

      End Do

      write (3,*) Null

      write (3,*)
      write (3,*)  " prtcl_ax: Particle No, Axial force "
      write (3,*)

      write (6,*)
      write (6,*)  " prtcl_ax: Particle No, Axial force "
      write (6,*)

      write (3,*) Nprtcl

      Do i=1,Nprtcl
        write (6,102) i,forcex(i)
        write (3,102) i,forcex(i)
      End Do

      write (3,*) Null
      write (3,*) 

c-------------------------------------
c END OF SOLVING THE INTEGRAL EQUATION
c-------------------------------------

      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to draw streamlines"
      write (6,*) " 0 to quit "
      write (6,*) " ----------"
      Read  (5,*) menu

      If(menu.eq.0) Go to 99

  92  Continue

      write (6,*) " Enter the maximum number of points along"
      write (6,*) " a streamline before inquiring "
      write (6,*) " ------------------------------"
      read  (5,*) Mstr

      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 2 to select the second-order RK method"
      write (6,*) " 4 to select the fourth-order RK method"
      write (6,*) " 0 to quit"
      write (6,*) " -------- "
      read  (5,*) IRK

      If(IRK.eq.0) Go to 99

      write (6,*) " Enter the spatial step"
      write (6,*) " ----------------------"
      read  (5,*) spat_step

      open (9,file="strml.dat")

  22  Continue

      read (9,*) X00,Y00

      If(abs(X00-99).lt.0.0000001) Go to 99   ! no more streamlines
      If(abs(Y00-99).lt.0.0000001) Go to 99

      call prtcl_ax_str
     +
     +   (Iflow
     +   ,X00,Y00
     +   ,Mstr
     +   ,IRK
     +   ,spat_step
     +   ,N_strml   ! number of points along the streamline
     +   )

      if(N_strml.eq.0) Go to 22  ! compute another streamline

c---
c Print out the streamline
c and its image
c---

      write (6,*) " prtcl_ax: one streamline with ",N_strml," points "
      write (6,*) "           completed"

      write (7,104) N_strml
      Do I=1,N_strml
        write (7,104) I,xstr(I),ystr(I)
      End Do

      write (7,104) N_strml 
      Do I=1,N_strml
        write (7,104) I,xstr(I),-ystr(I)
      End Do

c-------------
      Go to 22   ! compute another streamline
c-------------

c-----
c Done
c-----

  99  Continue

      write (3,*) Null
      write (7,*) Null

c---
c record parameters
c---

      write (3,169) Iflow
      write (3,170) visc
      write (3,171) NGL
      write (3,180) Uprtcl
      write (3,181) Uinf
      write (3,182) a11,a12
      write (3,183) a21,a22
      write (3,184) Max1,Max2
      write (3,186) Iprec
      write (3,986) Ireduce
      write (3,188) Nblocks
      write (3,189) (Lump(i),i=1,Nblocks)
      write (3,190) Niter
      write (3,990) tol
      write (3,192) Nprtcl

      Do i=1,Nprtcl
       write (3,193) 
     + j,xcntr(i),ycntr(i),axis1(i),axis2(i),NE(i),Itp(i)
      End Do

c---
c close files
c---

      close (3)
      close (7)
      close (9)

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
 180  Format (1x,"Uprtcl    = ",f10.5)
 181  Format (1x,"Uinf      = ",f10.5)
 182  Format (1x,"a11,a12   = ",f10.5,1x,f10.5)
 183  Format (1x,"a21,a22   = ",f10.5,1x,f10.5)
 184  Format (1x,"Max1, Max2= ",i2,1x,i2)
 186  Format (1x,"Iprec     = ",i1)
 986  Format (1x,"Ireduce   = ",i1)
 188  Format (1x,"Nblocks   = ",i3)
 189  Format (1x,"lump      = ",50(1x,i2))
 190  Format (1x,"Niter     = ",i5)
 990  Format (1x,"Tol       = ",f15.10)
 192  Format (1x,"Nprtcl    = ",i5)
 193  Format (1x,i2,1x,5(1x,f10.5),1x,2(1x,i2))

 888  Format (1x,"Iter: ",i3," Point corr: ",f15.10)

      Stop
      End
