      program body_ax

c==========================================
c FDLIB, CFDLAB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c==========================================

c---------------------------------------------------
c Axisymmetric streaming (uniform)
c potential flow past a stationary body.
c
c This program solves
c an integral equation of the second kind
c for the disturbance harmonic potential
c over the body contour
c and computes streamlines
c
c The contour of the body in the xy upper half-plane
c consists of a Number of SeGments (NSG)
c
c The segments may be straight lines or circular arcs.
c
c Each segment is discretized into a number of elements
c
c
c Symbols:
c --------
c
c NSG: Number of segments defining the body contour
c
c NE(i): Number of elements on the ith segment 
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
c T0(i): angle subtended from the center of a circular element
c        at ith collocation point
c
c arel(i):  axisymmetric surface area of ith element
c
c phi: disturbance potential at collocation points
c
c dphidn0: normal derivative of disturbance potential 
c                 at collocation points 
c
c dphids0: tangential derivative of potential 
c                 at collocation points 
c
c cp:	pressure coefficient at collocation points
c
c NGL: Number of Gaussian points for integration 
c      over each element
c
c Icross: Index for stopping the computation of the streamlines
c         Default value is 0
c         If Icross = 1, computation stops when a streamline
c         crosses the yz plane
c
c Notes:
c ------
c
c Normal vector points into the flow
c
c----------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension    NE(10),   RT(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,129),yw(10,129),tw(10,129)
      Dimension           phi(10,128)

      Dimension      X0(1280),  Y0(1280),  T0(1280),s0(1280)
      Dimension dphidn0(1280),velt(1280),veln(1280),cp(1280)
      Dimension    arel(1280)
      Dimension    tnX0(1280),tnY0(1280)
      Dimension    vnX0(1280),vnY0(1280)

      Dimension AL(1280,1280),BL(1280),SOL(1280)  ! for the linear system

      Dimension ZZ(20),WW(20)

      Dimension xstr(900),ystr(900) ! for streamlines

c--------------
c common blocks
c--------------

      common/xxx01/Iflow,NSG,NGL,NE,Itp
      common/xxx02/xw,yw,tw
      common/xxx03/actis,xcntr,ycntr

      common/xxx04/Vx,cr,Xlvr,Ylvr

      common/xxx05/X0,Y0,T0,S0,dphidn0
      common/xxx06/tnx0,tny0,vnx0,vny0,arel

      common/xxx07/xcenter,ycenter

      common/xxx08/xwmin,ywmin,xwmax,ywmax

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi  = 3.14159 265358  D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

      Null = 0
      None = 1
      Ntwo = 2

      zero  = 0.0D0
      five  = 5.0D0
      fivem =-5.0D0

c------
c input
c------

  94  Continue 

      write (6,*) 
      write (6,*) " Enter: "
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 50 for flow past a sphere"
      write (6,*) " 51 for flow past a triangular torus"
      write (6,*) "-------------------------------------"

      read (5,*) Iflow

      If(Iflow.eq.0) Go to 99

      If(Iflow.ne.50.and.Iflow.ne.51) then
        write (6,*) "Invalid selection; please try again"
        Go to 94
      End If

c-----------
c initialize
c-----------

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

      call body_ax_geo (Ncl)

c---------------------
c open streamline file
c---------------------

      open (9,file="body_ax.str")

c-------------------------
c Printing session
c
c Print boundary geometry
c to accompany streamlines
c-------------------------

      Ntot = 0
      Do j=1,NSG
       Ntot = Ntot+NE(j)+1
      End Do
      write (9,102) Ntot

      Do j=1,NSG
c      write (9,102) NE(j)+1
c      write (6,102) NE(j)+1
       Do i=1,NE(j)+1
        write (9,102) i,XW(j,i),YW(j,i)
c       write (6,102) i,XW(j,i),YW(j,i)
       End Do
      End Do

c---
c print the point vortex
c---

      write (9,102) None
      write (9,102) None,Xlvr,Ylvr

c---------------------------
c Display collocation points
c and disturbance velocity
c---------------------------

      write (6,*)
      write (6,*) Ncl," Collocation points "

c     write (6,*)
c     write (6,*) "Collocation points"
c     write (6,*) "and disturbance normal velocity"
c     write (6,*)
c
c     Do i=1,Ncl
c       write (6,102) i,X0(i),Y0(i),dphidn0(i)
c     End Do

c--------
c Prepare 
c--------

      call gauss_legendre (NGL,ZZ,WW)

c---------------------------------------------
c Generate the linear system
c for the potential at the collocation points
c
c Generate the influence matrix
c consisting of integrals of the
c single-layer potential
c
c Compute the rhs by evaluating the dlp
c-------------------------------------------------------
	
      Do 31 i=1,Ncl    ! loop over collocation points

       BL(i) = 0.0D0     ! initialize the right-hand side

c      write (6,*)        " Collocation point :",i

       j = 0              ! counter

       Do k=1,NSG         ! loop over segments

        If(Itp(k).eq.2) then
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

         j = j+1
         Ising = 0
         If(i.eq.j) Ising = 1

         call body_ax_sdlp
     +
     +     (X0(i),Y0(i),T0(i)
     +     ,X1,Y1,T1
     +     ,X2,Y2,T2
     +     ,NGL
     +     ,Ising
     +     ,Itp(k)
     +     ,rad,xcnt,ycnt
     +     ,QQQ
     +     ,WWW
     +     )

         AL(i,j) = WWW

         BL(i) = BL(i) + QQQ*dphidn0(j)

         End Do

        End Do

        AL(i,i) = AL(i,i) - 0.5D0

  31  Continue

c------------------------
c Solve the linear system
c------------------------

      write (6,*) " body_ax: Size of the linear system:",Ncl

      Isym_g = 0   ! system is not symmetric
      Iwlpvt = 1   ! pivoting enabled

      call gel
     +
     +   (Ncl,AL,BL,SOL
     +   ,Isym_g
     +   ,Iwlpvt
     +   ,Deter
     +   ,Istop
     +   )

      write (6,*) " body_ax: Linear system solved"

c-------------------
c Display the system
c-------------------

c     Do I=1,Ncl
c      write (6,101) (AL(i,j),j=1,Ncl),BL(I),SOL(I)
c     End Do

c------------------------
c Distribute the solution
c------------------------

      k = 0        ! counter

      Do i=1,NSG
       Do j=1,NE(i)
        k = k+1
        phi(i,j) = SOL(k)
       End Do
      End Do

c------------------------------------------------
c Compute the tangential disturbance velocity
c             by taking tangential derivative of phi
c
c Compute the pressure coefficient cp
c             drag force
c------------------------------------------------

c     write (6,*)
c     write (6,*) "Tang and normal tot vel; press coeff"
c     write (6,*)

      Forcex = 0.0

      k = 0       ! counter

      Do i=1,NSG
        Do j=1,NE(i)

         k = k+1

c--------------------------------
c tangential disturbance velocity
c computed by differentiating phi
c--------------------------------

         If(j.eq.1) then           ! forward differencing

          dphids0 = (phi(i,j+1)-phi(i,j))/(s0(k+1)-s0(k)) 

         Else If(j.eq.NE(i)) then  ! backward differencing

          dphids0 = (phi(i,j)-phi(i,j-1))/(s0(k)-s0(k-1)) 

         Else                      ! second-order differencing

          g1 = phi(i,j-1)
          g2 = phi(i,j)
          g3 = phi(i,j+1)
          h1 = s0(k-1)
          h2 = s0(k)
          h3 = s0(k+1)
          aa = ((g3-g2)/(h3-h2)-(g1-g2)/(h1-h2))/(h3-h1)
          bb = (g3-g2)/(h3-h2)-aa*(h3-h2)
          dphids0 = bb

         End If

c---------------
c total velocity
c---------------

         velx = dphidn0(k)*vnx0(k)+dphids0*tnx0(k)
         vely = dphidn0(k)*vny0(k)+dphids0*tny0(k)

         velx = velx + Vx        ! add incident flow

c---
c add line vortex ring
c---

         Iopt = 1

         call lvr_fs
     +
     +      (Iopt
     +      ,X0(k),Y0(k)
     +      ,Xlvr,Ylvr
     +      ,ulvr,vlvr
     +      ,psi
     +      )

         velx = velx + cr*ulvr
         vely = vely + cr*vlvr

c-------------------------------
c tangential and normal velocity
c and pressure coefficient 
c-------------------------------

         velt(k) = velx*tnx0(k)+vely*tny0(k)
         veln(k) = velx*vnx0(k)+vely*vny0(k)

c------------
c axial force
c------------

         cp(k) = 1.0D0-velt(k)**2/Vx**2  

         Forcex = Forcex + cp(k)*vnx0(k)*arel(k)

c        write (6,102) k,velt(k),veln,cp(k)

        End Do
      End Do

      write (6,*) 
      write (6,888) Forcex

c-------------------
c print the solution
c-------------------

      open (3,file="body_ax.out")

c     write (6,*)
c     write (6,*) "x, y, arc length, phi,"
c     write (6,*) "tangential and normal velocity, cp"
c     write (6,*)

      k = 0

      Do i=1,NSG

c       write (6,102) NE(i)
        write (3,102) NE(i)

        Do j=1,NE(i)
         k = k+1
         write (3,102) j,X0(k),Y0(k),s0(k),phi(i,j)
     +                 ,velt(k),veln(k),cp(k)
c        write (6,102) j,X0(k),Y0(k),s0(k),phi(i,j)
c    +                 ,velt(k),veln(k),cp(k)
        End Do

      End Do

      write (3,102) Null

      close (3)

c-----------------------------------------------------------

c     write (6,*)
c     write (6,*) "          MENU"
c     write (6,*)
c     write (6,*) " Enter 0 to quit"
c     write (6,*) "       1 to draw streamlines"
c     write (6,*) " ---------------------------"
c     read  (5,*) menu

c     If(menu.eq.0) Go to 99

c     write (6,*)
c     write (6,*) " Enter maximum number of points before "
c     write (6,*) "                  stopping or inquiring"
c     write (6,*) " --------------------------------------"
c     read  (5,*) Mstr

c     write (6,*)
c     write (6,*) " Enter 0 to stop when number is exceeded"
c     write (6,*) "       1 to inquire for continuation"
c     write (6,*) " --------------------------------------"
c     read  (5,*) Isc

c     write (6,*) 
c     write (6,*) "Integration method"
c     write (6,*) 
c     write (6,*) " Enter 0 to quit"
c     write (6,*) "       2 for the second-order Runge-Kutta"
c     write (6,*) "       4 for the fourth-order Runge-Kutta"
c     write (6,*) " ----------------------------------------"
c     read  (5,*) IRK

c     If(IRK.eq.0) Go to 99

c     write (6,*) 
c     write (6,*) " Enter the spatial step"
c     write (6,*) " ----------------------"
c     read  (5,*) Dl
c
c-----------------------

      menu  = 1
      Isc   = 0
      Mstr  = 300
      Dl    = 0.05
      IRK   = 2

c--------------------------
c prepare for crossing test
c--------------------------

      If(Iflow.eq.50) then          ! sphere
        Icross = 1
      Else If(Iflow.eq.51) then     ! torus
        Icross = 0
      End If

c--------------------------
c begin drawing streamlines
c--------------------------

      write (6,*)
      read  (4,*)

  22  Continue

      read (4,*) X00,Y00

c------------------------------
c will stop if X00=99 or Y00=99
c------------------------------

      If(abs(X00-99).lt.0.0000001) then
        close (4)
        Go to 97
      End If

      If(abs(Y00-99).lt.0.0000001) then
        close (4)
        Go to 97
      End If

      Xcross = X00   ! to be used for crossing check

c------

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
     +   (phi
     +   ,DphiDn0
     +   ,X00,Y00
     +   ,Ux1,Uy1
     +   )

c     write (6,104) L,X00,Y00,Ux1,Uy1

      step = Dl/sqrt(Ux1**2+Uy1**2)     ! set the frozen-time step

      Xsv = X00  ! save
      Ysv = Y00  ! save

c----------------------
      If(IRK.eq.2) then
c----------------------

        steph = 0.5*step
        X00 = Xsv + step * Ux1
        Y00 = Ysv + step * Uy1

        call velocity 
     + 
     +   (phi
     +   ,DphiDn0
     +   ,X00,Y00
     +   ,Ux2,Uy2
     +   )

        X00 = Xsv + steph*(Ux1+Ux2)
        Y00 = Ysv + steph*(Uy1+Uy2)

c---------------------------
      Else If(IRK.eq.4) then
c---------------------------

        steph = 0.5*step
        step6 = step/6.0

        X00 = Xsv + steph * Ux1
        Y00 = Ysv + steph * Uy1

        call velocity
     + 
     +   (phi
     +   ,DphiDn0
     +   ,X00,Y00
     +   ,Ux2,Uy2
     +   )

        X00 = Xsv + steph * Ux2
        Y00 = Ysv + steph * Uy2

        call velocity
     + 
     +   (phi
     +   ,DphiDn0
     +   ,X00,Y00
     +   ,Ux3,Uy3
     +   )

        X00 = Xsv + step * Ux3
        Y00 = Ysv + step * Uy3

        call velocity 
     + 
     +   (phi
     +   ,DphiDn0
     +   ,X00,Y00
     +   ,Ux4,Uy4
     +   )

        X00 = Xsv + step/6.0D0 * (Ux1+2.0*Ux2+2.0*Ux3+Ux4)
        Y00 = Ysv + step/6.0D0 * (Uy1+2.0*Uy2+2.0*Uy3+Uy4)

c-----------
      End If
c-----------

      K = K+1
      L = L+1

c---------------------------
c test for x=0 plane crossing
c---------------------------

      If(Icross.eq.1) then
       test = Xcross*X00
       If(test.lt.0) then
        write (6,*) " Crossed the x=0 plane: I will stop"
        Go to 21
       End If
      End If

c-------------------------
c test for sphere crossing
c-------------------------

      If(Iflow.eq.50) then

      crosss = Dsqrt((X00-xcenter)**2
     +             +(Y00-ycenter)**2)

      If(crosss.lt.actis(1)) then
        L = L-1
        Go to 23
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
      write (6,*) "Continue this streamline ?"
      write (6,*) 
      write (6,*) "Enter 0 for no, 1 for yes"
      write (6,*) "---------------------------"
      read  (5,*) Icon

      If(Icon.eq.1) Go to 20
     
c--------------------
c End of a streamline
c--------------------
  
  21  Continue

      xstr(L) = X00
      ystr(L) = Y00

  23  Continue

c-------------------------
c Print out the streamline
c-------------------------

      write (6,*) " One streamline with ",L," points completed"

      write (9,104) L  ! number of points

      Do I=1,L
        write (9,104) I,xstr(I),ystr(I)
      End Do

      Go to 22

c--------------------------
c End of streamline drawing
c--------------------------

  97  Continue

      write (9,104) Null
      close (9)

c-----
c Done
c-----

  99  Continue

 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
 150  Format (1X,10(1X,f10.5))
 888  Format (1X," Axial Force: ",f10.5)

      Stop
      End
