      program flow_1d

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------------
c Steady unidirectional flow through a tube
c with arbitrary cross-section
c
c This program solves
c an integral equation for the boundary distribution
c of the shear stress corresponding to
c the homogeneous component of the velocity field,
c and then evaluates the velocity at a specified point
c in the tube.
c
c The flow occurs along the z axis.
c
c The boundary of the tube in the xy plane
c consists of a Number of SeGments (NSG)

c The segments can be straight lines or circular arcs.
c
c Each segment is discretized into a number of elements
c
c
c Symbols:
c --------
c
c NSG:     Number of segments
c
c NE(i):   Number of elements on ith segment 
c
c RT(i):   Stretch ratio of elements on the ith segment 
c
c elml(i): length of the ith element
c
c Itp(i):  Index for shape of the ith segment:
c          1 for a straight segment
c          2 for a circular arc
c
c Ncl:     number of collocation points
c
c Xw, Yw:  end-nodes of elements on all segments
c
c X0, Y0:  coordinates of collocation points
c
c t0:      angle measured around the center a circular element
c
c dudn:    normal derivative of velocity at segment elements

c U0:      disturbance velocity at collocation points
c
c NGL:     Number of Gaussian points for integration 
c          over each element
c
c pg:      negative of the pressure gradient
c
c
c Notes:
c ------
c
c The normal vector points into the fluid
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(10),Itp(10),Actis(10)
      Dimension xcntr(10),ycntr(10)

      Dimension xw(10,300),yw(10,300),tw(10,300)
      Dimension          dudn(10,300)

      Dimension   X0(900),  Y0(900),T0(900),S0(0:900)
      Dimension   U0(900)
      Dimension tnX0(900),tnY0(900)
      Dimension vnx0(900),vny0(900)
      Dimension elml(900)

      Dimension AL(500,500),BL(500),SOL(500)  ! for the linear system

      Dimension ZZ(20),WW(20)    ! for the quadrature

c--------------
c common blocks
c--------------

      common/VEL00/Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw,U0
      common/VEL03/actis,xcntr,ycntr

      common/rectangle/sizex,sizey

      common/xxx04/visc,pg,orgx,orgy
      common/xxx05/X0,Y0,T0,S0
      common/xxx06/tnx0,tny0,vnx0,vny0
      common/xxx07/elml
      common/xxx08/xwmin,ywmin,xwmax,ywmax
      common/xxx10/xcnt,ycnt

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi = 3.1415926535897932384D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

      Null = 0
      None = 1
      Ntwo = 2

      zero  = 0.0D0
      five  = 5.0D0
      fivem =-5.0D0

c------------------------------------------

  93  Continue

      write (6,*) 
      write (6,*) " SELECT THE TUBE GEOMETRY"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) "  1 for a circular tube"
      write (6,*) "  2 for an elliptical tube"
      write (6,*) "  3 for a rectangular tube"
      write (6,*) "  4 for a triangular tube"
      write (6,*) " 10 to read data from file: contour.dat"
      write (6,*) "  0 to quit"
      write (6,*) " ----------"

      read (5,*) Iflow

      if(Iflow.eq.0) Go to 99

      if(Iflow.ne.1.and.Iflow.ne.2
     +             .and.Iflow.ne.3
     +             .and.Iflow.ne.4
     +             .and.Iflow.ne.10
     +  ) then
        write (6,*)
        write (6,*) " flow_1d: this option is not available"
        write (6,*) "          Please try again"
        write (6,*)
        Go to 93
      end if

c---------------------------------------------------------

c------------------------------
c boundary-element distribution
c------------------------------

       call flow_1d_geo
     +
     +   (Ncl   ! number of collocation points
     +   )

c-------------------------
c printing session
c
c print the boundary elements
c-------------------------

      open (9,file="flow_1d.elm",status='unknown')

      write (6,*)
      write (6,*) " Boundary geometry: element end-nodes"
      write (6,*) " ------------------------------------"
      write (6,*)

      Do j=1,NSG

       write (9,102) NE(j)+1
       write (6,102) NE(j)+1

       Do i=1,NE(j)+1
        write (9,102) i,XW(j,i),YW(j,i)
        write (6,102) i,XW(j,i),YW(j,i)
       End Do

      End Do

c---
c Display collocation points
c and disturbance velocity
c---

      write (6,*)
      write (6,*) "flow_1d:  ",Ncl," collocation points"
      write (6,*)
      write (6,*) "Coordinates and disturbance velocity:"
      write (6,*)

      Do i=1,Ncl
        write (6,102) i,X0(i),Y0(i),U0(i)
      End Do

c---------
c  prepare
c---------

      call gauss_leg (NGL,ZZ,WW)

c--------------------------
c Generate the linear system
c for the normal derivative
c of the velocity corresponding to the
c particular solution
c at the collocation points
c--------------------------

c--------------------------
c Generate the influence matrix
c consisting of integrals
c of the single-layer potential
c
c Compute the rhs by evaluating
c the double-layer potential
c--------------------------

      Do i=1,Ncl    ! loop over collocation points

c      write (6,*)        " Collocation point :",i

       BL(i) = -0.5D0*U0(i)     ! right-hand side

c      WWWsum = 0.0D0       ! for testing
c      QQQsum = 0.0D0       ! for testing

       j = 0              ! counter

       Do k=1,NSG         ! loop over segments

        if(Itp(k).eq.2) then
         rad  = actis(k)
         xcnt = xcntr(k)
         ycnt = ycntr(k)
        end if

        Do l=1,NE(k)       ! loop over elements

         X1 = XW(k,l)
         Y1 = YW(k,l)
         T1 = TW(k,l)

         X2 = XW(k,l+1)
         Y2 = YW(k,l+1)
         T2 = TW(k,l+1)

         j = j+1
         Ising = 0

         if(i.eq.j) Ising = 1

         call flow_1d_sdlp
     +
     +     (Iflow
     +     ,X0(i),Y0(i),t0(i)
     +     ,X1,Y1,T1
     +     ,X2,Y2,T2
     +     ,NGL
     +     ,Ising
     +     ,Itp(k)
     +     ,rad,xcnt,ycnt
     +     ,QQQ
     +     ,WWW
     +     )

         AL(i,j) = QQQ 
         BL(i)   = BL(i) + WWW*U0(j)   ! right-hand side

c        QQQsum = QQQsum + QQQ   ! for testing
c        WWWsum = WWWsum + WWW   ! for testing

         End Do

        End Do

c       write (6,100) i,QQQsum,WWWsum   ! for testing

      End Do

c------------------------
c Solve the linear system
c------------------------

      write (6,*)
      write (6,*) " flow_1d: size of the linear system: ",Ncl
      write (6,*)

      Isym_g = 0   ! system is not symmetric
      Iwlpvt = 1   ! pivoting enabled

      call gel
     +
     +    (Ncl
     +    ,AL,BL,SOL
     +    ,Isym_g
     +    ,Iwlpvt
     +    ,Det
     +    ,Istop
     +    )

c---
c display the linear system
c---
c
c      Do I=1,Ncl
c        write (6,101) (AL(I,J),J=1,Ncl),BL(I),SOL(I)
c        write (6,101) BL(i)
c      End Do

c------------------------
c distribute the solution
c------------------------

      Ic = 0        ! counter

      Do i=1,NSG

       Do j=1,NE(i)
        Ic = Ic+1
        dudn(i,j) = SOL(Ic)
       End Do

      End Do

c--------------------------------------------------
c Compute: 
c
c  1) the wall shear stress
c  2) the boundary integral of the boundary traction
c  3) the flow rate
c
c  Followed by a printing session
c--------------------------------------------------

      open (3,file="flow_1d.trc",status='unknown')

      write (6,*)
      write (6,*)  "     x, y, arc length, shear stress"
      write (6,*)

      Ic   = 0          ! collocation point counter
      bit  = 0.0D0      ! boundary integral of the traction
      flrt = 0.0D0

      fc = 0.50D0*pg/visc

      Do i=1,NSG

        write (6,102) NE(i)
        write (3,102) NE(i)

        Do j=1,NE(i)

         Ic = Ic+1

         dudnt = dudn(i,j)                    ! total normal derivative
     +     -fc*((x0(Ic)-orgx)*vnx0(Ic)
     +         +(y0(Ic)-orgy)*vny0(Ic) )

         shear_stress = visc*dudnt

         bit = bit + shear_stress*elml(Ic)

         flrtt = dudnt*U0(Ic)*visc/pg
     +         +pg*((x0(Ic)-orgx)**3*vnx0(Ic)
     +             +(y0(Ic)-orgy)**3*vny0(Ic))/(12.0D0*visc)

         flrt  = flrt + flrtt*elml(Ic)

         write (6,102) j,X0(Ic),Y0(Ic),S0(Ic),shear_stress
         write (3,102) j,X0(Ic),Y0(Ic),S0(Ic),shear_stress

        End Do

      End Do

c---
c normalize and print
c---

      bit  = bit/pi
      flrt = 8.0D0*flrt/pi

      write (3,102) Null
      write (3,*)
      write (3,106) bit
      write (3,*)
      write (3,107) flrt
      write (3,*)

      write (6,*)
      write (6,106) bit
      write (6,*)
      write (6,107) flrt

c-------------------------------

      open (8,file="flow_1d.vel")
      open (4,file="flow_1d.prof")

 92   Continue

      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 1 to evaluate the velocity at a point"
      write (6,*) " 2 for a 3D velocity graph"
      write (6,*) " 0 to continue"
      write (6,*) " -------------"
      read  (5,*) menu

c-----------------------
      if(menu.eq.1) then
c-----------------------

c Evaluate the velocity at a point

       write (6,*)
       write (6,*) " Please enter the x and y coordinates "
       write (6,*) " of the velocity point"
       write (6,*) " ---------------------"

       read  (5,*) Xvel,Yvel

       call velocity
     +   
     +     (Xvel,Yvel
     +     ,dudn
     +     ,vel)

       vel = vel                            ! total velocity
     +     - 0.25D0*pg/visc*((Xvel-orgx)**2
     +                      +(Yvel-orgy)**2
     +                      )

       write (6,108) vel
       write (8,101) Xvel,Yvel,vel

      Go to 92

c-----------------------
      else if(menu.eq.2) then
c-----------------------

c     If(Iflow.eq.3)  then ! rectangular tube

      write (4,*) Ncl
      Do i=1,Ncl
       write (4,110) X0(i),Y0(i)
      End Do

      Ndivx = 32
      Ndivy = 32

      sizex = 0.6D0
      sizey = 0.6D0

      Dx = 2.0*sizex/Ndivx
      Dy = 2.0*sizey/Ndivy

      write (4,*) Ndivx,Dx,sizex
      write (4,*) Ndivy,Dy,sizey

      Do j=1,Ndivy+1

       Yvel = (j-1.0)*Dy-sizey

       Do i=1,Ndivx+1

       Xvel = (i-1.0)*Dx-sizex

       If(i.eq.1.or.i.eq.Ndivx+1.or
     +   .j.eq.1.or.j.eq.Ndivy+1) then

        vel = 0.0D0

       else

       call velocity
     +   
     +     (Xvel,Yvel
     +     ,dudn
     +     ,vel)

       vel = vel                            ! total velocity
     +     - 0.25D0*pg/visc*((Xvel-orgx)**2
     +                      +(Yvel-orgy)**2
     +                      )
 
        if(vel.lt.0) vel = 0.0D0

        end if

c       write (6,120) i,j,Xvel,Yvel,vel
        write (4,110) vel

       End Do
      End Do

c     End If

c-----------
      End If
c-----------

      Go to 93

c-----
c done
c-----

   99 Continue

      write (3,104) Null
      write (8,104) Null
      write (9,104) Null

      close (8)
      close (9)
      close (3)

 100  Format (1x,i3,20(1x,f15.10))
 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
 106  Format (" Integral of the boundary traction : ",f15.10)
 107  Format (" Flow rate : ",f15.10)
 108  Format (" Velocity: ",f10.5)
 110  Format (4(1x,f9.5))
 120  Format (1x,i3,1x,i3,4(1x,f9.5))
 150  Format (1X,10(1X,f10.5))

      Stop
      End
