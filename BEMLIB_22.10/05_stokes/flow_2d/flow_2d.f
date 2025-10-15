      program flow_2d

c-----------------------------------------
c FDLIB - BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------
c Two-dimensional Stokes flow in a
c domain with arbitrary geometry
c
c This program solves the integral equation
c of the first kind
c for the traction around the boundaries,
c and then computes streamlines
c originating from specified points
c in the fluid.
c
c The boundary of the flow in the xy plane,
c or one period of it,
c consists of a number of segments
c that can be straight lines or circular arcs.
c
c Each segment is discretized into boundary
c elements with corresponding shapes
c
c SYMBOLS:
c --------
c
c NSG: Number of segments
c
c NE(i): Number of elements on ith segment 
c
c RT(i): Stretch ratio of elements on the ith segment 
c
c Itp(i): Index for the shape of the ith segment:
c         1 for a straight segment
c         2 for a circular arc
c
c (Xe, Ye):  end-nodes of elements
c (Xm, Ym):  mid-nodes of elements
c
c (X0, Y0):  coordinates of collocation points
c
c t0:	     angle subtended from a circular segment center
c
c (ux0, uy0): perturbation velocity components at collocation points
c
c fx, fy, Cartesian components of the traction
c
c NGL: Number of Gaussian points for integration over each element
c
c shrt:	shear rate of incident flow at the wall
c
c delta: pressure gradient of incident flow
c
c ftn: tangential component of the traction
c
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)
      Dimension fx(10,200),fy(10,200)

      Dimension   X0(900),  Y0(900),T0(900),S0(900)
      Dimension tnX0(900),tnY0(900)
      Dimension  ux0(900), uy0(900)

      Dimension elml(1280)

      Dimension ftn(10,200)             ! shear stress

      Dimension AL(900,900),BL(900),SOL(900)  ! for the linear system

      Dimension eig(900),eigt(900)      ! left and right eigenvectors of AL

      Dimension ZZ(20),WW(20)           ! Gauss-Legendre quadrature

      Dimension xstr(800),ystr(800)     ! for streamlines

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw,fx,fy,ux0,uy0
      common/xxx03/actis,xcntr,ycntr
      common/xxx04/visc,shrt,delta
      common/xxx05/X0,Y0,T0,S0
      common/xxx06/tnx0,tny0,vnx0,vny0
      common/xxx07/elml
      common/xxx08/xwmin,ywmin,xwmax,ywmax
      common/xxx10/xcnt,ycnt

      common/flow_91/RL,Uslip

      common/wwww/wall,rotation

      common/ZZWW/ZZ,WW

      common/ppii/pi,pih,pi2,pi4,pi6,pi8

c----------
c constants
c----------

      pi  = 3.14159 265358 D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi
      pi6 = 6.0D0*pi
      pi8 = 8.0D0*pi

      Null = 0
      None = 1

      zero = 0.0D0

c------
c input
c------

  94  Continue

      write (6,*) 
      write (6,*) " Please enter:"
      write (6,*) 
      write (6,*) "  1 for a circular cavity on a plane wall"
      write (6,*) " 11 for a rectangular cavity on a plane wall"
      write (6,*) " 41 for a circular protrusion on a plane wall"
      write (6,*) " 51 for a rectangular protrusion on a plane wall"
      write (6,*) " 71 for a rectangle above a plane wall"
      write (6,*) " 81 for a circular cylinder above a plane wall"
      write (6,*) " 91 for shear flow past a periodic array"
      write (6,*) "    circular cylinders"
      write (6,*) " 92 for shear flow past a periodic array"
      write (6,*) "    of rectangular cylinders"
      write (6,*) 
      write (6,*) "  0 to quit"
      write (6,*) "-----------"

      read  (5,*) Iflow

c------
c traps
c------

      If(Iflow.eq.0) Go to 99

      if    (Iflow.ne.1 
     +  .and.Iflow.ne.11
     +  .and.Iflow.ne.41
     +  .and.Iflow.ne.51
     +  .and.Iflow.ne.71
     +  .and.Iflow.ne.81
     +  .and.Iflow.ne.91
     +  .and.Iflow.ne.92
     +  ) then

        write (6,*)
        write (6,*) " flow_2d: this selection is not avaliable;"
        write (6,*) "          please try again"

        Go to 94

      End  If

c-----------
c initialize
c-----------

c----------------------
c  Compute:
c
c  the element distribution
c  the collocation points
c  the disturbance velocity 
c      at the collocation points
c
c  upper wall is located at y = wall
c----------------------

      call flow_2d_geo
     +
     +  (Ncl
     +  ,Iwall
     +  ,Itry
     +  )

      write (6,*)
      write (6,*) " flow_2d: ",Ncl," collocation points"
      write (6,*)

c---------------------
c open streamline file
c---------------------

      open (9,file="flow_2d.str")

c----------------------------------------
c print elements to accompany streamlines
c----------------------------------------

c---
c print the lower-wall nodes
c
c NSGL: number of boundary segments
c---

      if(Iflow.eq.1)  NSGL = 4
      if(Iflow.eq.11) NSGL = 5
      if(Iflow.eq.41) NSGL = 4
      if(Iflow.eq.51) NSGL = 5
      if(Iflow.eq.71) NSGL = 4         ! square
      if(Iflow.eq.81) NSGL = 1         ! circle
      if(Iflow.eq.91) NSGL = 1
      if(Iflow.eq.92) NSGL = 4

      Ntot = 0
      Do j=1,NSGL
       Ntot = Ntot+NE(j)+1
      End Do
      write (9,102) Ntot

      Do j=1,NSGL
c      write (9,102) NE(j)+1
c      write (6,102) NE(j)+1
       Do i=1,NE(j)+1
        write (9,102) i,XW(j,i),YW(j,i)
c       write (6,102) i,XW(j,i),YW(j,i)
       End Do
      End Do

c---
c print the upper wall nodes
c---

      If(Iwall.eq.1) then

       Ntot = 0
       Do j=NSGL+1,NSG
        Ntot = Ntot+NE(j)+1
       End Do
       write (9,102) Ntot

       Do j=NSGL+1,NSG
c       write (9,102) NE(j)+1
c       write (6,102) NE(j)+1
        Do i=1,NE(j)+1
         write (9,102) i,XW(j,i),YW(j,i)
c        write (6,102) i,XW(j,i),YW(j,i)
        End Do
       End Do

      End If

c-----------------------------------
c Display the collocation points
c and print the disturbance velocity
c-----------------------------------

c      write (6,*)
c      write (6,*) " flow_2d: Coll points and dist velocity"
c      write (6,*)

c      Do i=1,Ncl
c       write (6,102) i,X0(i),Y0(i),ux0(i),uy0(i)
c      End Do

c---------------
c prepare to run
c---------------

      call gauss_leg (NGL,ZZ,WW)

      Ncl2 = 2*Ncl

      cf_slp = -1.0D0/(pi4*visc)
      cf_dlp =  1.0D0/ pi4

c-------------------------------------------
c Generate the linear system
c for the traction at the collocation points
c
c The influence matrix consists 
c of integrals of the single-layer potential
c
c Compute the rhs by evaluating the
c double-layer potential
c------------------------------------------

      Do 31 I=1,Ncl    ! loop over collocation points

c      write (6,*)    " flow_2d: collocating at point :",I

       NI = Ncl+I

c---
c right-hand side
c note that: u_0 = - u_inf
c---

       If(   Iflow.eq.71
     +   .or.Iflow.eq.81
     +   .or.Iflow.eq.91
     +   .or.Iflow.eq.92
     +   ) then
         BL(I)  = ux0(I)
         BL(NI) = uy0(I)
       Else
         BL(I)  = 0.5D0*ux0(I)
         BL(NI) = 0.5D0*uy0(I)
       End If

c---
c run over elements
c---

       J = 0              ! element counter

       Do k=1,NSG         ! loop over segments

        If(Itp(k).eq.2) then   ! for arcs
         rad  = actis(k)
         xcnt = xcntr(k)
         ycnt = ycntr(k)
        End If

        Do L=1,NE(k)       ! loop over elements

         X1 = XW(K,L)
         Y1 = YW(K,L)
         T1 = TW(K,L)

         X2 = XW(K,L+1)
         Y2 = YW(K,L+1)
         T2 = TW(K,L+1)

         J = J+1

         Ising = 0
         If(I.eq.J) Ising = 1   ! subtract off the singularity

         call flow_2d_sdlp
     +
     +      (Iflow
     +      ,X0(I),Y0(I),t0(I)
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,Qxx,Qxy
     +      ,Qyx,Qyy
     +      ,Wxx,Wyx
     +      ,Wxy,Wyy
     +      )

         NJ = Ncl +J

         AL(I, J) = cf_slp*Qxx     ! x component of the BIE
         AL(I,NJ) = cf_slp*Qyx     ! x component of the BIE

         AL(NI, J) = cf_slp*Qxy     ! y component of the BIE
         AL(NI,NJ) = cf_slp*Qyy     ! y component of the BIE

         If(   Iflow.eq.71
     +     .or.Iflow.eq.81
     +     .or.Iflow.eq.91
     +     .or.Iflow.eq.92
     +   ) then
          BL(I)  = BL(I)
          BL(NI) = BL(NI)
         Else
          BL(I)  = BL(I) -(ux0(J)*Wxx+uy0(J)*Wyx)*cf_dlp
          BL(NI) = BL(NI)-(ux0(J)*Wxy+uy0(J)*Wyy)*cf_dlp
         End If

         End Do      ! run over elements
        End Do       ! run over segments

  31  Continue      ! run over collocation points

c----------------------------------------
c In the case of shear flow over an array of cylinders
c (Iflow = 91, 92)
c add one more unknown corresponding to the
c slip velocity
c
c and introduce one more equation determining
c the force exerted on each cylinder:
c
c the integral of the traction over a cylinder
c should be equal to visc*shrt*RL
c
c----------------------------------------

      if(Iflow.eq.91.or.Iflow.eq.92) then
      
       Ncl1  = Ncl +1
       Ncl21 = Ncl2+1

c---
c one more column numbered Ncl2+1
c---

       Do i=1,Ncl
        AL(i,Ncl21) = 1.0D0
       End Do

       Do i=Ncl1,Ncl21
        AL(i,Ncl21) = 0.0D0
       End Do

c---
c one more equation
c---

       Do i=1,Ncl
        AL(Ncl21,i) = elml(i)
c       write (6,*) i,elml(i)
       End Do

       Do j=Ncl1,Ncl21
        AL(Ncl21,j) = 0.0D0
       End Do

       BL(Ncl21) = visc*shrt*RL
c      write (6,*) BL(Ncl21)

       write (6,*) " flow_2d: extended MLS constructed"

      End If

c----------------------------------------
c Flow past a block
c
c Compute the eigenvectors
c of the influence matrix
c and regularize the system
c----------------------------------------

      if(Iflow.eq.71.or.Iflow.eq.81) then

      Ic = 0

      If(Iflow.eq.71) then    ! flow past a square

       Do i=1,NSG
        Do j=1,NE(i)
         Ic = Ic+1
         IcNcl = Ic+Ncl
         dxx = xw(I,J+1)-xw(I,J)
         dyy = yw(I,J+1)-yw(I,J)
         dss = Dsqrt(dxx**2+dyy**2)
         eigt(Ic)    = -dyy
         eigt(IcNcl) =  dxx
         eig (Ic)    = -dyy/dss
         eig (IcNcl) =  dxx/dss
        End Do
       End Do

      End If

      If(Iflow.eq.81) then   ! flow past a cylinder

       Do j=1,NE(1)
          Ic  = Ic+1
          IcNcl = Ic+Ncl
          elal = actis(1)*(tw(1,j+1)-tw(1,j))
          eigt(Ic)    = Dcos(T0(Ic))
          eigt(IcNcl) = Dsin(T0(Ic))
          eig (Ic)    =    eigt(Ic)/elal
          eig (IcNcl) = eigt(IcNcl)/elal
       End Do

      End If

c---
c normalize the eigenvectors
c---

      sum  = 0.0D0
      sumt = 0.0D0

      Do i=1,Ncl2
       sum  = sum + eig(i)**2
       sumt = sumt+eigt(i)**2
      End Do
      sum  = Dsqrt(sum)
      sumt = Dsqrt(sumt)
      Do i=1,Ncl2
        eig(i) =  eig(i)/sum
       eigt(i) = eigt(i)/sumt
      End Do

c---
c check identities
c---

      Go to 333

      Do I=1,Ncl2
       sum1 = 0.0D0
       sum2 = 0.0D0
       sum3 = 0.0D0
       Do J=1,Ncl2
        sum1 = sum1+ eig(J)*AL(I,J)
        sum2 = sum1+eigt(J)*AL(J,I)
        sum3 = sum3+eigt(J)*BL(J)
       End Do
       write (6,105) I,sum1,sum2,sum3
      End Do

 333  Continue

c---
c orthogonal projection
c---

      sum = 0.0D0
      Do i=1,Ncl2
        sum = sum+eigt(i)*BL(i)
      End Do
      Do i=1,Ncl2
        BL(i) = BL(i)-sum*eigt(i)
      End Do

      Do i=1,Ncl2
       sum = 0.0D0
       Do j=1,Ncl2
        sum = sum+eigt(j)*AL(j,i)
       End Do
       Do j=1,Ncl2
         AL(i,j) = AL(i,j)-sum*eigt(j)
       End Do
      End Do

c-----------
      End If    ! end of regularization
c-----------

c----------------------------------
c Set the size of the linear system
c----------------------------------

      Nsys = Ncl2

      if(Iflow.eq.91.or.Iflow.eq.92) then
        Nsys = Ncl2+1      ! the slip velocity is an additional unknown
      end if
      
      write (6,*)
      write (6,*) " flow_2d: Size of the linear system: ",Nsys

c------------------------
c Solve the linear system
c------------------------

      Isym_g = 0   ! system is not symmetric
      Iwlpvt = 1   ! pivoting enabled

      call gel
     +
     +    (Nsys
     +    ,AL,BL,SOL
     +    ,Isym_g
     +    ,Iwlpvt
     +    ,Deter
     +    ,Istop
     +    )

       If(Istop.eq.1) stop

c------------------------------------
c uncomment to view the linear system
c or the solution
c------------------------------------
 
c      write (6,*) 
c      write (6,*)  "flow_2d: linear system"
c      write (6,*) 

c      Do I=1,Nsys
c        write (6,101) (AL(i,j),J=1,Nsys),BL(i),SOL(i)
c      End Do
 
c      write (6,*) 
c      write (6,*)  "flow_2d: solution vector"
c      write (6,*) 

c      Do i=1,Nsys
c       write (6,101) SOL(i)
c      End Do
 
c-------------------------
c  Distribute the solution
c-------------------------

      K = 0        ! element counter

      Do I=1,NSG
       Do J=1,NE(i)

        K = K+1
        fx(I,J) = SOL(K)
        fy(I,J) = SOL(K+Ncl)

       End Do
      End Do

c----------------------------------------
c For shear flow over an array of cylinders
c (Iflow = 91, 92)
c recover the slip velocity and compute
c the drift velocity
c----------------------------------------

      If(Iflow.eq.91.or.Iflow.eq.92) then

        Uslip = SOL(Nsys)

        Ic = 0

        sum = 0.0D0

        Do i=1,NSG
         Do j=1,NE(i)
          Ic = Ic+1
          sum  = sum + Y0(Ic)*SOL(Ic)*elml(Ic)
         End Do
        End Do

        Udrift = Uslip + sum/RL

        write (6,*) 
        write (6,134) Uslip,Udrift
        write (6,*) 

      End If

c----------------------------------
c Compute the tangential components
c of the traction over the element
c----------------------------------

      K = 0        ! element counter

      Do I=1,NSG
       Do J=1,NE(I)

        K = K+1
        ftn(I,J) = fx(I,J)*tnx0(K)+fy(I,J)*tny0(K)

       End Do
      End Do

c---
c Add shear stress of incident flow
c over the wall
c if desired
c---
c
c     const = visc*shrt
c
c     Do J=1,NE(NSG)
c       Fty(NSG,J) = const + Fty(NSG,J) 
c       Ftz(NSG,J) = const + Ftz(NSG,J)
c     End Do

c-------------------
c print the solution
c-------------------

c     write (6,*)
c     write (6,*)  " x, y, arc_length, fx, fy, shear stress "
c     write (6,*)
c     if(Iflow.ne.91.and.Iflow.ne.92
c    +              .and.Iflow.ne.71.and.Iflow.ne.81) then
c       write (6,*)  " (disturbance variables)"
c     end if
c     write (6,*)  " -----------------------------------"

      open (3,file="flow_2d.out")

      Ic = 0          ! element counter

      Do k=1,NSG

c       write (6,107) NE(i)
        write (3,107) NE(i)

        Do j=1,NE(k)
         Ic = Ic+1
c        write (6,107) j,X0(Ic),Y0(Ic),S0(Ic)
c    +                  ,fx(k,j),fy(k,j),ftn(k,j)
         write (3,107) j,X0(Ic),Y0(Ic),S0(Ic)
     +                  ,fx(k,j),fy(k,j),ftn(k,j)
        End Do

      end Do

      write (3,107) Null

      close (3)

c----------------------------------------------
c Flow over a circular protrusion
c
c Compute force and torque
c
c TORK is torque with respect to the origin
c TORQ is torque with respect to the center 
c----------------------------------------------

      if(Iflow.eq.41) then

      DRAG = 0.0D0
      TORK = 0.0D0
      TORQ = 0.0D0

      Ic = NE(1)   ! collocation point counter

      Do k=2,3
       Do j=1,NE(k)
        Ic = Ic+1
        elal = actis(k)*(tw(k,j+1)-tw(k,j))
        elal = abs(elal)
        DRAG = DRAG + elal * fx(k,j)
        TORK = TORK + elal * (X0(Ic)*fy(k,j)-Y0(Ic)*fx(k,j) )
        TORQ = TORQ + elal *((X0(Ic)-xcnt)*fy(k,j)
     +                      -(Y0(Ic)-ycnt)*fx(k,j))
       End Do
      End Do

      write (6,150) DRAG,TORK,TORQ
      write (6,108) DRAG,TORK,TORQ

      End If

c-------------------------
c Flow past a rectangular protrusion
c
c Compute force and torque
c-------------------------

      If(Iflow.eq.51) then

c     open (2,file="flow_2d.out1")

      DRAG = 0.0D0
      TORK = 0.0D0
      TORQ = 0.0D0

      Ic = 0  ! collocation point counter

      Do k=2,4
       Do j=1,NE(k)
        Ic = Ic+1
        elal = sqrt((xw(k,j+1)-xw(k,j))**2   ! element arc length
     +             +(yw(k,j+1)-yw(k,j))**2)
        fcon = elal * fx(k,j)
        tcon = elal * (X0(Ic)*fy(k,j)-Y0(Ic)*fx(k,j) )
        qcon = elal *((X0(Ic)-xcnt)*fy(k,j)
     +               -(Y0(Ic)-ycnt)*fx(k,j))
        DRAG = DRAG + fcon
        TORK = TORK + tcon
        TORQ = TORQ + qcon
c       write (6,107) J,X0(Ic),Y0(Ic),S0(Ic),fcon,tcon,qcon
c       write (2,107) J,X0(Ic),Y0(Ic),S0(Ic),fcon,tcon,qcon
       End Do
      End Do

      write (6,150) DRAG,TORK,TORQ
      write (6,108) DRAG,TORK,TORQ
c     close (2)

      End If

c-------------------------
c Flow past a square
c
c Compute force and torque
c-------------------------

      If(Iflow.eq.71) then

c     open (2,file="flow_2d.out1")

      DRAG = 0.0D0
      TORK = 0.0D0
      TORQ = 0.0D0

      Ic = 0  ! collocation point counter

      Do k=1,4
       Do j=1,NE(k)
        Ic = Ic+1
        elal = sqrt((xw(k,j+1)-xw(k,j))**2   ! element arc length
     +             +(yw(k,j+1)-yw(k,j))**2)
        fcon = elal * fx(k,j)
        tcon = elal * (X0(Ic)*fy(k,j)-Y0(Ic)*fx(k,j) )
        qcon = elal *((X0(Ic)-xcnt)*fy(k,j)
     +               -(Y0(Ic)-ycnt)*fx(k,j))
        DRAG = DRAG + fcon
        TORK = TORK + tcon
        TORQ = TORQ + qcon
c       write (6,107) J,X0(Ic),Y0(Ic),S0(Ic),fcon,tcon,qcon
c       write (2,107) J,X0(Ic),Y0(Ic),S0(Ic),fcon,tcon,qcon
       End Do
      End Do

      write (6,150) DRAG,TORK,TORQ
      write (6,108) DRAG,TORK,TORQ
c     close (2)

      End If

c------------------------------
c Flow past a circular cylinder
c
c Compute force and torque
c------------------------------

      If(Iflow.eq.81) then

      DRAG = 0.0D0
      TORK = 0.0D0
      TORQ = 0.0D0

      Ic = 0        ! collocation point counter

      Do j=1,NE(1)
       Ic = Ic+1
       elal = actis(1)*(tw(1,j+1)-tw(1,j))
       elal = abs(elal)
       DRAG = DRAG + elal * fx(1,j)
       TORK = TORK + elal * (X0(Ic)*fy(1,j)-Y0(Ic)*fx(1,j) )
       TORQ = TORQ + elal *((X0(Ic)-xcnt)*fy(1,j)
     +                     -(Y0(Ic)-ycnt)*fx(1,j))
      End Do

      write (6,150) DRAG,TORK,TORQ
      write (6,108) DRAG,TORK,TORQ

      End If

c-------------------------------

      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 to draw streamlines"
      write (6,*) " 2 to draw a velocity profile"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c-----------------------
      if(menu.eq.1) then    ! will draw streamlines
c-----------------------

c     write (6,*)
c     write (6,*) " The streamlines will be recorded in"
c     write (6,*) " file: flow_2d.str"
c     write (6,*)
c     write (6,*) " Enter the maximum number of points"
c     write (6,*) "       along a streamline before pausing"
c     write (6,*) " -------------------------------------- "
c     read  (5,*) Mstr

c     write (6,*)
c     write (6,*) " Enter 0 to stop when number is exceeded"
c     write (6,*) "       1 to inquire for continuation"
c     write (6,*) " --------------------------------------"
c     read  (5,*) Isc

c     write (6,*)
c     write (6,*) " Select the integration method"
c     write (6,*)
c     write (6,*) " Enter 0 to quit"
c     write (6,*) "       1 for the  first-order Runge-Kutta method"
c     write (6,*) "       2 for the second-order Runge-Kutta method"
c     write (6,*) "       4 for the fourth-order Runge-Kutta method"
c     write (6,*) " -----------------------------------------------"
c     read  (5,*) IRK

c     If(IRK.eq.0) Go to 99

c     write (6,*)
c     write (6,*) " Enter the spatial step"
c     write (6,*) " ----------------------"
c     read  (5,*) Dl
c
c-----------------------

      Mstr  = 600
      Isc   = 0
      IRK   = 2
      Dl    = 0.01

c--------------------------
c begin drawing streamlines
c--------------------------

      write (6,*)

      Icross = 1

      if(Iflow.eq.71) Icross = 0

      read  (4,*)  ! file opened in subroutine: flow_2d_geo

c------
c read the starting point of the streamline
c------

  22  Continue     ! will return here to draw another streamline

      read (4,*) X00,Y00                 ! 99 to quit

      If(abs(X00-99).lt.0.0000001) then
        close (4)
        Go to 97
      End If
      If(abs(Y00-99).lt.0.0000001) then
        close (4)
        Go to 97
      End If

      Xcross = X00  ! for crossing check
      Xclosd = X00  ! for closed streamline check
      Yclosd = Y00  ! for closed streamline check

c---

      L = 1     ! local counter for inquiry
      K = 1     ! total counter

  20  Continue

      xstr(L) = X00
      ystr(L) = Y00

c---------------
c integrate ODEs
c---------------

      call velocity 
     +
     +   (X00,Y00
     +   ,Ux1,Uy1
     +   )

c     write (6,104) L,X00,Y00,Ux1,Uy1

      step = Dl/sqrt(Ux1*Ux1+Uy1*Uy1)        ! time step

c---------------
c First-order Rk
c---------------

      if(IRK.eq.1) then

        X00 = X00 + step*Ux1
        Y00 = Y00 + step*Uy1

c----------------
c Second-order Rk
c----------------

      elseif(IRK.eq.2) then

        Xsv = X00      ! save
        Ysv = Y00

        steph = 0.5*step

        X00 = Xsv + step * Ux1
        Y00 = Ysv + step * Uy1

        call velocity 
     +
     +     (X00,Y00
     +     ,Ux2,Uy2
     +     )

        X00 = Xsv + steph*(Ux1+Ux2)
        Y00 = Ysv + steph*(Uy1+Uy2)

c----------------
c Fourth-order RK
c----------------

      elseif(IRK.eq.4) then

        Xsv = X00      ! save
        Ysv = Y00

        X00 = Xsv + 0.5*step * Ux1
        Y00 = Ysv + 0.5*step * Uy1

        call velocity
     +
     +    (X00,Y00
     +    ,Ux2,Uy2
     +    )

        X00 = Xsv + 0.5*step * Ux2
        Y00 = Ysv + 0.5*step * Uy2

        call velocity
     +
     +     (X00,Y00
     +     ,Ux3,Uy3
     +     )

        X00 = Xsv + step * Ux3
        Y00 = Ysv + step * Uy3

        call velocity 
     +
     +     (X00,Y00
     +     ,Ux4,Uy4
     +     )

        X00 = Xsv + step * (Ux1+2.0*Ux2+2.0*Ux3+Ux4)/6.0
        Y00 = Ysv + step * (Uy1+2.0*Uy2+2.0*Uy3+Uy4)/6.0

c-----------
      End If
c-----------

      K = K+1
      L = L+1

c------------------
c crossing test
c remove if desired
c------------------

      If(Icross.eq.1) then
       test = Xcross*X00
       If(test.lt.0) then
c       write (6,*) " flow_2d: Crossed the yz plane"
        Go to 21
       End If
      End If

c------------------------
c closed streamline check
c------------------------

      test = sqrt((X00-Xclosd)**2+(Y00-Yclosd)**2)

      tol = Dl/2.0

      If(test.lt.tol) then
        X00 = Xclosd
        Y00 = Yclosd
        write (6,*) " flow_2d: Streamline closed"
        Go to 21
      End If

c-----------------------
c window crossing checks
c-----------------------

      If(X00.lt.xwmin) Go to 21
      If(X00.gt.xwmax) Go to 21
      If(Y00.lt.ywmin) Go to 21
      If(Y00.gt.ywmax) Go to 21

c---------------------
c point capacity check
c---------------------

      If(K.lt.Mstr) Go to 20    ! continue the streamline

      If(Isc.eq.0)  Go to 21    ! abandon this streamline

      K = 1    ! reset local counter

      write (6,*) 
      write (6,*) "Continue this streamline ?"
      write (6,*) 
      write (6,*) "Enter 0 for NO, 1 for YES"
      write (6,*) "-------------------------"
      read  (5,*) Icon

      If(Icon.eq.1) Go to 20

c--------------------
c End of a streamline
c--------------------

  21  Continue

      xstr(L) = X00
      ystr(L) = Y00

c-------------------------
c Print out the streamline
c-------------------------

      write (6,*) " One streamline with ",L," points completed"

      If(L.gt.3) then

       write (9,104) L      ! number of points

       Do I=1,L
        write (9,104) I,xstr(I),ystr(I)
       End Do

      End If

      Go to 22   ! Go back for another streamline

c--------------------------
c End of streamline drawing
c--------------------------

  97  Continue

      write (9,104) Null
      close (9)

c------------------------------------------
      End If   ! end of drawing streamlines
c------------------------------------------

      If(menu.eq.2) then

c-----------------------------
c Will draw a velocity profile
c-----------------------------

      open (7,file="flow_2d.prof")

      write (6,*)
      write (6,*) " Will draw a velocity profile"
      write (6,*) "      and compute the flow rate"
      write (6,*)
      write (6,*) " Please enter Xmin Xmax"
      write (6,*) " ----------------------"
      read  (5,*) Xmin,Xmax

      write (6,*) " Please enter Ymin Ymax"
      write (6,*) " ----------------------"
      read  (5,*) Ymin,Ymax

      write (6,*) " How many points for the profile ?"
      write (6,*) " ---------------------------------"
      read  (5,*) Nprof

      Nprof1 = Nprof+1
      write (6,104) Nprof1

c---
c steps
c---

      DX0 = (Xmax-Xmin)/(Nprof1-1.0)
      DY0 = (Ymax-Ymin)/(Nprof1-1.0)

c---
c normal vector
c---

      vnx =  DY0
      vny = -DX0

c---
c profiling
c---

      X00 = Xmin
      Y00 = Ymin

      flow = 0.0

      Do i=1,Nprof1

        call velocity 
     +
     +    (X00,Y00
     +    ,Ux,Uy
     +    )

        write (6,104) I,X00,Y00,Ux,Uy
        write (7,104) I,X00,Y00,Ux,Uy

c---
c flow rate by the trapezoidal rule
c---

        Q = Ux*vnx+Uy*vny
        If(i.gt.1) flow = flow + 0.5*(Q+Qsave)
        Qsave = Q

        X00 = X00 + DX0
        Y00 = Y00 + DY0

      End Do

      write (6,*)
      write (6,*) " Flow rate = ",flow
      write (6,*)

      close (7)

c------------------------------------------
      End If     ! and of drawing a profile
c------------------------------------------

      Go to 94

c-----
c Done
c-----

  99  Continue

 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
 105  Format (1x,i3,9(1x,f15.10))
 107  Format (1x,i3,20(1x,f11.5))
 108  Format (20(1x,f12.8))
 134  Format (" Slip  velocity = ",f15.10,/
     +       ," Drift velocity = ",f15.10
     +       )
 150  Format ("DRAG:   ",F15.10,/
     +       ,"TORQUE: ",F15.10,/
     +       ,"TORK:   ",F15.10
     +       )

      Stop
      End
