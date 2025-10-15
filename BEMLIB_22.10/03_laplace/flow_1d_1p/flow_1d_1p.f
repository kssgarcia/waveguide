      program flow_1d_1p

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c Steady unidirectional shear flow over a
c periodic array of cylinders.
c
c This program solves
c an integral equation for the boundary distribution
c of the shear stress,
c and then evaluates the velocity at a specified point
c in the flow.
c
c The flow occurs along the z axis, parallel to the
c cylinder genrators.
c
c The boundary of the cylinderd in the xy plane
c consists of a Number of SeGments (NSG)

c The segments can be straight lines or circular arcs.
c
c Each segment is discretized into a number of elements
c
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
c elml(i): length of the ith element
c
c Itp(i): Index for shape of the ith segment:
c         1 for a straight segment
c         2 for a circular arc
c
c Xe, Ye:  end-nodes of elements on a segment
c Xm, Ym:  mid-nodes of elements on a segment
c Xw, Yw:  end-nodes of elements on all segments
c
c X0, Y0:  coordinates of collocation points
c
c t0: angle measured around the center a circular element
c
c dudn: normal derivative of velocity at segment elements

c u0: disturbance velocity at collocation points
c
c NGL: Number of Gaussian points for integration 
c      over each element
c
c
c Notes:
c ------
c
c The normal vector points into the fluid
c
c
c Capacity:
c ---------
c
c  10 Segments,
c 128 elemens over each segment
c
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(10),   RT(10),Itp(10),actis(10)
      Dimension xcntr(10),ycntr(10)

      Dimension xw(10,129),yw(10,129),tw(10,129)
      Dimension          dudn(10,128)

      Dimension   x0(1280),  y0(1280),t0(1280)
      Dimension   s0(1280),  u0(1280)
      Dimension tnX0(1280),tnY0(1280)
      Dimension vnX0(1280),vnY0(1280)
      Dimension elml(1280)

c---
c for the linear system
c---

      Dimension AL(1280,1280),BL(1280),SOL(1280)

c---
c various
c---

      Dimension ZZ(20),WW(20)

c---
c common blocks
c---

      common/VEL00/RL,visc,Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw
      common/VEL03/actis,xcntr,ycntr
      common/VEL04/x0,y0,t0,s0,tnx0,tny0,vnx0,vny0,elml
      common/VEL05/u0
      common/VEL06/dudn

      common/REC01/recx,recy

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      Null = 0
      None = 1
      Ntwo = 2

      pi  = 3.14159265358D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi
      pih = 0.5D0*pi

c------
c input
c------

  93  Continue

      write (6,*) 
      write (6,*) " SELECT THE GEOMETRY"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) "  0 to quit"
      write (6,*) "  1 for circular cylinders"
      write (6,*) "  2 for elliptical cylinders"
      write (6,*) "  3 for rectangular cylinders"
      write (6,*) "  4 for triangular cylinders"
      write (6,*) " 10 for rectangular protrusion"
      write (6,*) " -----------------------------"

      read (5,*) Iflow

      If(Iflow.eq.0) Go to 99

      If(Iflow.ne.1.and.Iflow.ne.2
     +             .and.Iflow.ne.3
     +             .and.Iflow.ne.4
     +             .and.Iflow.ne.10
     +  ) then
        write (6,*)
        write (6,*) " This is option is not available"
        write (6,*) " Please try again"
        write (6,*)
        Go to 93
      End If

c---------------------------------------------------------

c------------------------------
c Boundary element distribution
c------------------------------

      call flow_1d_1p_geo
     +
     +   (shrt
     +   ,Ncl
     +   )

c-------------------------
c Printing session
c
c print boundary elements
c-------------------------

      open (9,file="flow_1d_1p.elm")

      Go to 334

      write (6,*)
      write (6,*) " flow_1d_1p: Element end-nodes:"
      write (6,*) " -----------------------------"
      write (6,*)

      Do j=1,NSG

       write (9,102) NE(j)+1
       write (6,102) NE(j)+1

       Do i=1,NE(j)+1
        write (9,102) i,XW(j,i),YW(j,i)
        write (6,102) i,XW(j,i),YW(j,i)
       End Do

      End Do

 334  Continue

c---------------------------
c Display collocation points
c and undisturbed velocity
c---------------------------

      Go to 333

      write (6,*)
      write (6,*) "flow_1d_1p:  ",Ncl," collocation points"
      write (6,*)
      write (6,*) "flow_1d_1p: collocation point coordinates"
      write (6,*) "            and disturbance velocity"
      write (6,*) "------------------------------------------"

      Do i=1,Ncl
        write (6,102) i,X0(i),Y0(i),u0(i)
      End Do

 333  Continue

c----------------------------------------
c  PREPARE TO SOLVE THE INTEGRAL EQUATION
c----------------------------------------

      call gauss_legendre (NGL,ZZ,WW)

c--------------------------
c  Generate the linear system
c  for the normal derivative
c  of the velocity 
c  at the collocation points
c--------------------------

c--------------------------
c Generate the influence matrix
c consisting of integrals of the single-layer potential
c--------------------------

      Do 31 ic=1,Ncl      ! loop over collocation points

c      write (6,*)        " Collocation point :",i

       BL(ic) = visc*u0(ic)      ! right-hand side

       jc = 0              ! element counter

       Do k=1,NSG         ! loop over segments

        If(Itp(k).eq.2) then
         rad  = actis(k)
         xcnt = xcntr(k)
         ycnt = ycntr(k)
        End If

        Do L=1,NE(k)       ! loop over segment-elements

         X1 = XW(k,L)
         Y1 = YW(k,L)
         T1 = TW(k,L)

         X2 = XW(k,L+1)
         Y2 = YW(k,L+1)
         T2 = TW(k,L+1)

         jc = jc+1

         Ising = 0

         If(ic.eq.jc) Ising = 1

         call flow_1d_1p_slp
     +
     +      (Iflow
     +      ,RL
     +      ,X0(ic),Y0(ic),t0(ic)
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,rad,xcnt,ycnt
     +      ,QQQ
     +      )

         AL(ic,jc) = QQQ 

         End Do

        End Do

  31  Continue

c----------------------------------------
c Add one more unknown corresponding to the
c slip velocity
c
c Introduce one more equation determining
c the force exerted on each cylinder:
c
c the integral of the traction over a cylinder
c should be equal to visc*shrt*RL
c
c----------------------------------------

       Ncl1 = Ncl+1

c---
c one more column
c---

       Do i=1,Ncl
        AL(i,Ncl1) = -visc
       End Do

c---
c one more equation
c---

       Do j=1,Ncl
        AL(Ncl1,j) = elml(j)  ! element arc length
       End Do

       BL(Ncl1) = visc*shrt*RL

c------------------------
c SOLVE THE LINEAR SYSTEM
c------------------------

      Nsys = Ncl1    ! size of the linear system

      write (6,*)
      write (6,*) " flow_1p_1d: linear system size: ",Nsys
      write (6,*)

      Isym_g = 0   ! system is not symmetric
      Iwlpvt = 1   ! pivoting enabled

      call gel
     +
     +   (Nsys
     +   ,AL,BL,SOL
     +   ,Isym_g
     +   ,Iwlpvt
     +   ,Determ
     +   ,Istop
     +   )

c---
c Display the linear system
c---
c
c      Do I=1,Nsys
c        write (6,101) (AL(I,j),j=1,Nsys),BL(I),SOL(I)
c        write (6,101) BL(i)
c      End Do

c------------------------
c Distribute the solution
c------------------------

      K = 0        ! counter

      Do i=1,NSG
       Do j=1,NE(i)
        K = K+1
        dudn(i,j) = SOL(K)
       End Do
      End Do

c----------------------------------------
c recover the slip velocity and compute
c the drift velocity
c----------------------------------------

      Uslip = SOL(Nsys)

      sum = 0.0D0

      k = 0        ! element counter

      Do i=1,NSG
       Do j=1,NE(I)
        k = k+1
        sum = sum + Y0(k)*dudn(i,j)*elml(k)
       End Do
      End Do

      Udrift = Uslip+sum/RL

c--------------------------------------------------
c Compute and print:
c
c   shear stress
c   boundary integral of the total boundary traction
c--------------------------------------------------

      open (3,file="flow_1d_1p.trc",status='unknown')

      write (6,*)
      write (6,*)  "     x, y, arc length, shear stress"
      write (6,*)

      k = 0            ! counter
      bit = 0.0D0      ! boundary integral of the traction

      Do i=1,NSG

c       write (6,102) NE(i)
        write (3,102) NE(i)

        Do j=1,NE(i)
         k = k+1
         sstr = visc*dudn(i,j)
         bit = bit + sstr*elml(k)
c        write (6,102) j,x0(k),y0(k),s0(k),sstr,elml(k)
         write (3,102) j,x0(k),y0(k),s0(k),sstr
        End Do

      End Do

c-----------------
c printing session
c-----------------

      If(Iflow.eq.3) then
       a = 2.0D0*recx/RL
       b = 2.0D0*recy/RL
       Unorm   = - shrt/pi2 * log(Dcos(pih*(1.0D0-a)))
       Uslipn  = Uslip /Unorm
       Udriftn = Udrift/Unorm
       write (3,150) a,b,Uslip,Udrift,Uslipn,Udriftn,bit
       write (6,150) a,b,Uslip,Udrift,Uslipn,Udriftn,bit
      Else

      write (6,*)
      write (6,134) Uslip,Udrift
      write (6,*)
      write (6,106) bit

      write (3,102) Null
      write (3,*)
      write (3,106) bit
      write (3,*)

      End If

c-------------------------------

      open (8,file="flow_1d_1p.vel")

      write (6,*)
      write (6,*) "              MENU"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 to evaluate the velocity "
      write (6,*) "         at a point in the flow"
      write (6,*) "       2 to prepare a velocity matrix"
      write (6,*) " ------------------------------------"
      read  (5,*) menu

c     menu = 0

      If(menu.eq.0) Go to 93

c-----
c Evaluate the velocity at a point
c-----

      If(menu.eq.1) then

 92   Continue

       write (6,*)
       write (6,*) " Please enter the x and y coordinates "
       write (6,*) "        of the velocity point"
       write (6,*) "        99 to quit"
       write (6,*) " ----------------------------"

       read  (5,*) Xvel,Yvel

       If(Xvel.eq.99) Go to 99
       If(Yvel.eq.99) Go to 99

       call velocity 
     +
     +    (Xvel,Yvel
     +    ,shrt
     +    ,Uslip
     +    ,vel)

       vel = vel + Uslip              ! total velocity

       write (6,107) vel
       write (8,101) Xvel,Yvel,vel

      Go to 92

      close (8)

c-----
c Evaluate the velocity matrix
c-----

c---------------------------
      Else If(menu.eq.2) then
c---------------------------

      open (8,file="flow_1d_1p.vel")

      RLH = 0.5D0*RL

      Nx = 64
      Ny = 64
      Dx = 2.0*RL/Nx
      Dy = RL/Ny
c     ybottom = 0.001
      Dy = 2*RL/Ny
      ybottom = -RL-recy-0.01

      write (8,*) Nx
      write (8,*) Ny
      write (8,*) Dx
      write (8,*) Dy
      write (8,*) RL
      write (8,*) ybottom

      Do j=1,Ny+1
       Yvel = ybottom + (j-1.0D0)*Dy
       Do i=1,Nx+1
       Xvel = -RL + (i-1.0D0)*Dx
       call velocity 
     +
     +    (Xvel,Yvel
     +    ,shrt
     +    ,Uslip
     +    ,vel)
       write (8,*) vel
       End Do
      End Do

c---
c     Nx = 16
c     Ny = 16
c     Dx = 2.0*RL/Nx
c     Dy = RL/Ny
c     ybottom = -RL-recy-0.01
c     write (8,*) Nx
c     write (8,*) Ny
c     write (8,*) Dx
c     write (8,*) Dy
c     write (8,*) ybottom
c     Do j=1,Ny+1
c      Yvel = ybottom +(j-1.0D0)*Dy
c      Do i=1,Nx+1
c      Xvel = -RL + (i-1.0D0)*Dx
c      call velocity 
c    +
c    +    (Xvel,Yvel
c    +    ,shrt
c    +    ,Uslip
c    +    ,vel)
c      write (8,*) vel
c      End Do
c     End Do
c-----

      write (8,*) recx
      write (8,*) recy

      close (8)

c-----------
      End If
c-----------

c-----
c Done
c-----

   99 Continue

      write (3,104) Null
      write (8,104) Null
      write (9,104) Null

      close (9)
      close (3)

 100  Format (1x,i3,20(1x,f15.10))
 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
 106  Format (" Integral of the boundary traction : ",f15.12)
 107  Format (" Velocity: ",f10.5)
 134  Format (" Slip  velocity = ",f15.10,/
     +       ," Drift velocity = ",f15.10
     +       )
 150  Format (1X,10(2X,f12.8))

      Stop
      End
