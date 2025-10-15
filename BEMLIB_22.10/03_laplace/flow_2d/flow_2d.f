      program flow_2d

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c----------------------------------------------
c Two-dimensional potential flow
c in an arbitrary domain
c
c The program solves an integral equation
c for the harmonic potential
c and computes streamlines
c
c The boundary of the flow in the xy plane,
c or one period of the boundary,
c consists of a Number of SeGments (NSG)
c
c The segments can be straight lines or circular arcs.
c
c Each segment is discretized into a number of elements
c that conform with the segment geometry.
c
c
c Symbols:
c --------
c
c NSG: Number of segments
c
c NE(i): Number of elements on ith segment 
c
c RT(i): Stretch ratio of elements on ith segment 
c
c Itp(i): Index for shape of the ith segment:
c         1 for a straight segment
c         2 for a circular arc
c
c (Xe, Ye):  end-nodes of elements on a segment
c (Xm, Ym):  mid-nodes of elements on a segment
c (Xw, Yw):  end-nodes of elements on all segments
c
c (X0, Y0):  coordinates of collocation points
c
c  t0: angle subtended from the center of a circular element
c
c phi: disturbance potential at segment elemetns 
c
c dphidn0: normal derivative of potential 
c                 at collocation points 
c
c NGL: Number of Gaussian points for integration 
c      over each element
c
c Note:
c -----
c
c Normal vector points into the flow
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Xe(100),Ye(100),Te(100),Se(100)
      Dimension Xm(100),Ym(100),Tm(100),Sm(100)

      Dimension NE(10),RT(10),Itp(10),Actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)
      Dimension phi(10,200)

      Dimension X0(500),Y0 (500),T0(500),S0(0:500),dphidn0(500)
      Dimension tnX0(500),tnY0(500)
      Dimension vnX0(500),vnY0(500)

      Dimension AL(500,500),BL(500),SOL(500)  ! for the linear system

      Dimension ZZ(20),WW(20)

      Dimension xstr(900),ystr(900)    ! for streamlines

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw,phi,dphidn0
      common/xxx03/actis,xcntr,ycntr
      common/xxx04/Vx,Vy
      common/xxx05/X0,Y0,T0,S0
      common/xxx06/tnx0,tny0,vnx0,vny0
      common/xxx08/xwmin,ywmin,xwmax,ywmax

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

      Null = 0
      None = 1
      Ntwo = 2

      zero  =  0.0D0
      five  =  5.0D0
      fivem = -5.0D0

c------
c input
c------

  94  Continue

      write (6,*) 
      write (6,*) " Please enter:"
      write (6,*) 
      write (6,*) "  1 for a circular cavity on a plane wall "
      write (6,*) " 11 for a rectangular cavity on a plane wall"
      write (6,*) " 41 for a circular protrusion on a plane wall"
      write (6,*) "  0 to quit"
      write (6,*) " ----------"

      read (5,*) Iflow

c------
c traps
c------

      If(Iflow.eq.0) Go to 99

      If    (Iflow.ne.1.and.Iflow.ne.11
     +  .and.Iflow.ne.41
     +     ) then

        write (6,*)
        write (6,*) " This selection is not avaliable"
        write (6,*) " Please try again"
        Go to 94

      End  If

c-----------
c initialize
c-----------

      Itry  = 1      ! graphics

  93  Continue       ! graphics

c----------------------
c  Compute:
c
c  the element distribution
c  the collocation points
c  the normal derivative of the
c      disturbance potential
c      at the collocation points
c
c----------------------

      call flow_2d_geo (Ncl,Iwall,wall)

c---------------------
c open streamline file
c---------------------

      open (9,file="flow_2d.str")

c--------------------------------------------------
c Printing session
c
c Print boundary geometry to accompany streamlines
c--------------------------------------------------

c---
c print the lower-wall nodes
c---

      if(Iflow.eq.1)  NSGL = 4
      if(Iflow.eq.11) NSGL = 5
      if(Iflow.eq.41) NSGL = 4

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
c print the upper-wall nodes
c---

      If(Iwall.eq.1) then

      Ntot = 0

      Do j=NSGL+1,NSG
       Ntot = Ntot+NE(j)+1
      End Do

      write (9,102) Ntot

      Do j=NSGL+1,NSG
c      write (9,102) NE(j)+1
c      write (6,102) NE(j)+1
       Do i=1,NE(j)+1
        write (9,102) i,XW(j,i),YW(j,i)
c       write (6,102) i,XW(j,i),YW(j,i)
       End Do

      End Do

      End If

c---------------------------
c Display collocation points
c and disturbance normal velocity
c---------------------------

      write (6,*)
      write (6,*) "flow_2d: ",Ncl," Collocation points"

c     write (6,*) "Collocation points"
c     write (6,*) "and disturbance normal velocity"
c     write (6,*)
c
c     Do i=1,Ncl
c       write (6,102) i,X0(i),Y0(i),dphidn0(i)
c     End Do

c--------------------
c Read the quadrature
c--------------------

      call gauss_leg (NGL,ZZ,WW)

c------------------------------------------------------
c Generate the linear system
c for the potential at the collocation points
c
c Generate the influence matrix
c consisting of integrals of the single-layer potential
c
c Compute the rhs by evaluating the dlp
c------------------------------------------------------

      Do 31 I=1,Ncl    ! loop over collocation points

       BL(I) = 0.0D0     ! right-hand side

c      write (6,*)        " Collocation point :",I

       J = 0              ! counter

       Do k=1,NSG         ! loop over segments

        If(Itp(k).eq.2) then     ! arcs
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

         J   = J+1
         Ising = 0
         IF(I.eq.J) Ising = 1

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
     +      ,QQQ
     +      ,WWW
     +      )

         AL(I,J) = WWW

         BL(I) = BL(I)+QQQ*dphidn0(J)

         End Do

        End Do

        AL(I,I) = AL(I,I) - 0.5D0

  31  Continue

c------------------------
c solve the linear system
c------------------------

      write (6,*) " flow_2d: Size of the linear system:",Ncl

      Isym_g = 0   ! system is not symmetric
      Iwlpvt = 1   ! pivoting enabled

      call gel
     +
     +   (Ncl
     +   ,AL,BL,SOL
     +   ,Isym_g
     +   ,Iwlpvt
     +   ,Deter
     +   ,Istop
     +   )

      write (6,*) " flow_2d: Linear system solved"

c------------------------
c Print the linear system
c
c      Do I=1,Ncl
c        write (6,101) (AL(I,J),J=1,Ncl),BL(I),SOL(I)
c        write (6,101) BL(i)
c      End Do
c------------------------

c------------------------
c Distribute the solution
c------------------------

      k = 0        ! counter

      Do i=1,NSG
       Do J=1,NE(i)
        k = k+1
        phi(I,J) = SOL(k)
       End Do
      End Do

c------------------------------
c Record and print the solution
c------------------------------

      open (3,file="flow_2d.out")

c     write (6,*)
c     write (6,*) "x, y, arc length, dphi/dn, phi"
c     write (6,*)

      K = 0        ! counter

      Do i=1,NSG

c       write (6,102) NE(i)
        write (3,102) NE(i)

        Do J=1,NE(i)

         K = K+1
c        write (6,102) J,X0(K),Y0(K),S0(K),dphidn0(K),phi(i,j)
         write (3,102) J,X0(K),Y0(K),S0(K),dphidn0(K),phi(i,j)

        End Do

      End Do

      write (3,102) Null
      close (3)

c-------------------------------------------------------

c     write (6,*)
c     write (6,*) " Enter:"
c     write (6,*)
c     write (6,*) " 0 to quit"
c     write (6,*) " 1 to draw streamlines"
c     write (6,*) " 2 to draw a velocity profile"
c     write (6,*) " ----------------------------"
c     read  (5,*) menu

      menu = 1

      If(menu.eq.0) Go to 99

c-----------------------
      If(Menu.eq.1) then    ! will draw streamlines
c-----------------------

c     write (6,*)
c     write (6,*) " Enter maximum number of points before "
c     write (6,*) "                  stopping or inquiring"
c     write (6,*) " --------------------------------------"
c     read  (5,*) Mstr

c     write (6,*)
c     write (6,*) " Enter:"
c     write (6,*)
c     write (6,*) " 0 to stop when number is exceeded"
c     write (6,*) " 1 to inquire for continuation"
c     write (6,*) " -----------------------------"
c     read  (5,*) Isc

c     write (6,*)
c     write (6,*) " Enter:"
c     write (6,*)
c     write (6,*) " 0 to quit"
c     write (6,*) " 2 to select the second-order RK method"
c     write (6,*) " 4 to select the fourth-order RK method"
c     write (6,*) " ------------------------------------- "
c     read  (5,*) IRK

c     If(IRK.eq.0) Go to 99

c     write (6,*)
c     write (6,*) " Enter the spatial step"
c     write (6,*) " ----------------------"
c     read  (5,*) Dl
c
c-----------------------

      menu = 1
      Mstr = 300
      Isc  = 0
      IRK  = 2
      Dl   = 0.05D0

c--------------------------
c begin drawing streamlines
c--------------------------

      Icross = 1

      write (6,*)
      read  (4,*)

  22  Continue           ! will return here to draw another streamline

c---
c starting point
c---

      read (4,*) X00,Y00   ! will quit if X00=99 or Y00=99

      If(abs(X00-99).lt.0.0000001) then
        close (4)
        Go to 97
      End If
      If(abs(Y00-99).lt.0.0000001) then
        close (4)
        Go to 97
      End If

      Xcross = X00   ! for crossing check

c---

      L = 1          ! local counter for inquiry
      K = 1          ! total counter

  20  Continue

      xstr(L) = X00
      ystr(L) = Y00

c---------------
c integrate ODEs
c---------------

      call velocity 
     +
     +    (X00,Y00
     +    ,Ux1,Uy1
     +    )

c     write (6,104) L,X00,Y00,Ux1,Uy1

      step = Dl/Dsqrt(Ux1**2+Uy1**2)     ! set frozen-time step

c----------------------
      If(IRK.eq.2) then
c----------------------

        Xsv = X00  ! save
        Ysv = Y00  ! save

        steph = 0.5*step

        X00 = Xsv + step * Ux1
        Y00 = Ysv + step * Uy1

        call velocity 
     +
     +      (X00,Y00
     +      ,Ux2,Uy2
     +      )

        X00 = Xsv + steph*(Ux1+Ux2)
        Y00 = Ysv + steph*(Uy1+Uy2)

c----------------------
      Else If(IRK.eq.4) then
c----------------------

        Xsv = X00  ! save
        Ysv = Y00  ! save

        steph = 0.5*step
        step6 = step/6.0

        X00 = Xsv + steph * Ux1
        Y00 = Ysv + steph * Uy1

        call velocity
     +
     +     (X00,Y00
     +     ,Ux2,Uy2
     +     )

        X00 = Xsv + steph * Ux2
        Y00 = Ysv + steph * Uy2

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

        X00 = Xsv + step6 * (Ux1+2.0*Ux2+2.0*Ux3+Ux4)
        Y00 = Ysv + step6 * (Uy1+2.0*Uy2+2.0*Uy3+Uy4)

c-----------
      End If
c-----------

      K = K+1
      L = L+1

c---------------------------
c test for yz plane crossing
c---------------------------

      If(Icross.eq.1) then
       test = Xcross*X00
       If(test.lt.0) then
        write (6,*) " Crossed the yz plane: I will stop"
        Go to 21
       End If
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
      write (6,*) "  Continue this streamline ?"
      write (6,*) 
      write (6,*) "  Enter 0 for NO, 1 for YES"
      write (6,*) "---------------------------"
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

      write (9,104) L  ! number of points

      Do I=1,L
        write (9,104) I,xstr(I),ystr(I)
      End Do

      Go to 22   ! Go back for another streamline

c--------------------------
c End of streamline drawing
c--------------------------

  97  Continue

      write (9,104) Null
      close (9)

      Go to 94

c----------
      End If   ! end of drawing streamlines
c----------

      If(Menu.eq.2) then

      open (10,file="flow_2d.prof")

      write (6,*)
      write (6,*) " Will draw an x velocity profile"
      write (6,*) " across the y axis"
      write (6,*)

      write (6,*) " Please the x position"
      write (6,*) " ----------------------"
      read  (5,*)  Xprof

      write (6,*) " Please enter ymin ymax"
      write (6,*) " ----------------------"
      read  (5,*)  Ymin,Ymax

      write (6,*)
      write (6,*) " How many points for the profile ?"
      write (6,*) " ---------------------------------"
      read  (5,*) Nprof

      Nprof1 = Nprof+1

      Dy = (Ymax-Ymin)/(Nprof1-1.0D0)

      yprof = ymin

      write (6,104) Nprof1

      Do i=1,Nprof1

        call velocity 
     +
     +     (Xprof,Yprof
     +     ,Ux,Uy
     +     )

        Yprof = Yprof + Dy

        write (10,104) I,Xprof,Yprof,Ux,Uy
        write ( 6,104) I,Xprof,Yprof,Ux,Uy

      End Do

c----------
      End If     ! and of drawing a profile
c----------

  99  Continue

c-----
c Done
c-----

 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
 150  Format (1X,10(1X,f10.5))

      Stop
      End
