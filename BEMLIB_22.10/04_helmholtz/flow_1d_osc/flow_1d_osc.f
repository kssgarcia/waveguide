      program flow_1d_osc

c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c--------------------------------------------
c Oscillatory unidirectional flow through 
c straight tube with arbitrary cross-section
c or in the exterior of a rectilinear tube
c
c This program solves an integral equation
c for the boundary distribution of 
c the complex shear stress corresponding to
c the homogeneous component of the velocity,
c and then evaluates the velocity
c at a specified point in the flow
c
c The boundary of the tube in the xy plane
c consists of a Number of SeGments (NSG)

c The segments can be straight lines or circular arcs.
c
c Each segment is discretized into a number of elements
c
c
c Symbols:
c -------
c
c NSG:     Number of segments
c NE(i):   Number of elements on ith segment 
c RT(i):   Stretch ratio of elements on the ith segment 
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
c t0: angle subtended from the center a circular element
c
c dwdnr: real part of the normal derivative of the
c        amplitude of the velocity at segment elements
c
c dwdni: imaginary part of the normal derivative of the
c        amplitude of the velocity at segment elements
c
c w0: total or disturbance velocity amplitude
c      at collocation points
c
c NGL: Number of Gaussian points for integration 
c      over each element
c
c
c Notes:
c ------
c
c The normal vector must point toward
c the interior of the cylinder
c----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension NE(10),Itp(10)
      Dimension actis(10),xcntr(10),ycntr(10)

      Dimension xw(10,200),yw(10,200),tw(10,200)

      Dimension   X0(500),  Y0(500),T0(500),S0(0:500),w0(500)
      Dimension tnX0(500),tnY0(500)
      Dimension vnX0(500),vnY0(500)
      Dimension elml(500)

      Dimension dwdnr(10,200),dwdni(10,200)

c---
c for the linear system
c---

      Dimension AL(500,500),BL(500),SOL(500)

      Dimension ZZ(20),WW(20)

c--------------
c common blocks
c--------------

      common/VEL00/Iflow,NSG,NGL,NE,Itp
      common/VEL02/xw,yw,tw,w0
      common/VEL03/actis,xcntr,ycntr
      common/VEL05/dwdnr,dwdni

      common/xxx04/visc,den,omega
      common/xxx05/X0,Y0,T0,S0
      common/xxx06/tnx0,tny0,vnx0,vny0
      common/xxx07/elml
      common/xxx08/xwmin,ywmin,xwmax,ywmax
      common/xxx10/xcnt,ycnt

      common/ZZWW/ZZ,WW

c----------
c constants
c----------

      pi = 3.1415 92653 58979 32384 D0

      pih = 0.5D0*pi
      pi2 = 2.0D0*pi

      Null = 0
      None = 1
      Ntwo = 2

      zero  = 0.0D0
      five  = 5.0D0
      fivem =-5.0D0

c----------
c inquiries
c----------
       
 96   Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) "  1 for a  circular tube"
      write (6,*) "  2 for an elliptical tube"
      write (6,*) "  3 for a  rectangular tube"
      write (6,*) "  4 for a  triangular tube"
      write (6,*) "  0 to quit"
      write (6,*) "  ---------"

      read (5,*) Iflow

      If(Iflow.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for interior flow"
      write (6,*) " 2 for exterior flow"
      write (6,*) " 0 to quit"
      write (6,*) "----------"

      read (5,*) Inout

      If(Inout.eq.0) Go to 99

c--------
c prepare
c--------

      sign = 1.D0
      If(Inout.eq.2) sign = - sign

c----------------------
c element distribution
c----------------------

      call flow_1d_osc_geo (Ncl)

c-------------------------
c Printing session
c
c print element end-points 
c-------------------------

      open (9,file="flow_1d_osc.elm",status='unknown')

      write (6,*)
      write (6,*) " Element distribution:"
      write (6,*)

      Do j=1,NSG

       write (9,102) NE(j)+1
       write (6,102) NE(j)+1

       Do i=1,NE(j)+1
        write (9,102) i,XW(j,i),YW(j,i)
        write (6,102) i,XW(j,i),YW(j,i)
       End Do

      End Do

c----------------------------------
c Display collocation points
c and total or disturbance velocity
c----------------------------------

      write (6,*)
      write (6,*) " Collocation points and wall velocity:"
      write (6,*)

      Do i=1,Ncl
        write (6,102) i,X0(i),Y0(i),w0(i)
      End Do

c----------------------------------------
c  PREPARE TO SOLVE THE INTEGRAL EQUATION
c----------------------------------------

      call gauss_leg (NGL,ZZ,WW)

      delta = dsqrt(omega*den/visc)  ! Stokes layer thickness

      Ncl2 = 2*Ncl     ! size of the linear system

c----------------------------------------
c  Generate the linear system
c  for the real and imaginary part of the 
c  shear stress at the collocation points
c----------------------------------------

c------------------------------------------------------
c Generate the influence matrix of size 2*Ncl
c consisting of integrals of the real and
c imaginary part of the single-layer potential
c
c Compute the rhs by evaluating the dlp
c
c First set of Ncl equations will contain the real part
c           of the integral equation
c
c Second set of Ncl equations will contain the imaginary part
c           of the integral equation
c
c First set of Ncl unknowns will the real part of dwdn
c Second set of Ncl unknowns will the imaginary part of dwdn
c
c------------------------------------------------------

      Do 31 i=1,Ncl    ! loop over collocation points

c      write (6,*)        " Collocation point :",i

       incl = i+Ncl

       BL(i)    = - sign*0.5D0*w0(i)    ! right-hand side
       BL(incl) =        0.0D0          ! right-hand side

       If(menu.eq.3) BL(i) = -BL(i)

c      QQQrsum = 0.0D0      ! for testing
c      QQQisum = 0.0D0
c      WWWrsum = 0.0D0
c      WWWisum = 0.0D0

       J = 0              ! counter

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

         J = J+1
         Ising = 0

         If(I.eq.J) Ising = 1

         call flow_1d_osc_sdlp
     +
     +      (X0(I),Y0(I),t0(I)
     +      ,X1,Y1,T1
     +      ,X2,Y2,T2
     +      ,NGL
     +      ,Ising
     +      ,Itp(k)
     +      ,delta
     +      ,rad,xcnt,ycnt
     +      ,QQQr,QQQi      ! single layer
     +      ,WWWr,WWWi      ! double layer
     +      )

         jncl = j+Ncl

         AL(i,j)    =  QQQr  ! real part of the equation
         AL(i,jncl) = -QQQi  ! real part of the equation

         AL(incl,j)    =  QQQi  ! img part of the equation
         AL(incl,jncl) =  QQQr  ! img part of the equation

         BL(i)    = BL(i)   +WWWr*w0(J)   ! right-hand side
         BL(incl) = BL(incl)+WWWi*w0(J)   ! right-hand side

c        QQQrsum = QQQrsum + QQQr
c        QQQisum = QQQisum + QQQi
c        WWWrsum = WWWrsum + WWWr
c        WWWisum = WWWisum + WWWi

         End Do

        End Do

c       write (6,100) i,QQQrsum,WWWrsum,QQQisum,WWWisum

  31  Continue

c------------------------
c SOLVE THE LINEAR SYSTEM
c------------------------

      write (6,*) 
      write (6,*) "Size of the linear system:",Ncl2
      write (6,*) 

      Isym_g = 0        ! system is not symmetric
      Iwlpvt = 1        ! enables pivoting

      call gel
     +
     +    (Ncl2
     +    ,AL,BL,SOL
     +    ,Isym_g
     +    ,Iwlpvt
     +    ,Deter
     +    ,Istop
     +    )

c--------------------------
c display the linear system
c--------------------------

c      Do I=1,Ncl2
c        write (6,101) (AL(I,J),J=1,Ncl),BL(I),SOL(I)
c      End Do

c---------------------
c  DISTRIBUTE SOLUTION
c---------------------

      K = 0        ! counter

      Do i=1,NSG

       Do J=1,NE(i)
        K = K+1
        dwdnr(I,J) = SOL(K)          ! real part
        dwdni(I,J) = SOL(K+Ncl)      ! imaginary part
       End Do

      End Do

c--------------------------------
c Compute:
c         1) total boundary traction
c
c         2) boundary integral of
c            the total boundary traction
c--------------------------------

      open (3,file="flow_1d_osc.trc",status='unknown')

      write (6,*)
      write (6,*) " x, y, arcl, Re(shear stress), Im(shear stress)"
      write (6,*)

      k = 0       ! counter

      bitr = 0.0  ! real part of the boundary integral of traction
      biti = 0.0  ! imaginary part of the boundary integral of traction

      Do i=1,NSG

        write (6,102) NE(i)
        write (3,102) NE(i)

        Do j=1,NE(i)

         k = k+1

         shear_r = visc*dwdnr(i,j)                 ! shear stress
         shear_i = visc*dwdni(i,j)

         bitr = bitr + shear_r*elml(k)
         biti = biti + shear_i*elml(k)

         write (6,102) J,X0(K),Y0(K),S0(K),shear_r,shear_i
         write (3,102) J,X0(K),Y0(K),S0(K),shear_r,shear_i

        End Do

      End Do

c--------------------
c normalize and print
c--------------------

      bitr = visc*bitr/pi
      biti = visc*biti/pi

      write (3,102) Null
      write (3,*)
      write (3,106) bitr,biti
      write (3,*)

      write (6,*)
      write (6,106) bitr,biti

c-------------------------------

      open (8,file="flow_1d_osc.vel",status="unknown")

 92   Continue

      write (6,*)
      write (6,*) " Enter 0 to repeat"
      write (6,*) "       1 to evaluate the velocity "
      write (6,*) "         at a point in the flow"
      write (6,*) " ---------------------------------"
      read  (5,*) menu

      If(menu.eq.0) Go to 96

c-----------------------
      if(menu.eq.1) then    ! evaluate the velocity
c-----------------------

      write (6,*) 
      write (6,*) " Enter the x and y coordinates of the "
      write (6,*) " point where the velocity will be computed"
      write (6,*) " -------------------------------------"

      read  (5,*) Xvel,Yvel

      call velocity 
     +
     +   (Xvel,Yvel
     +   ,delta
     +   ,velr,veli
     +   )

      velr = sign*velr
      veli = sign*veli

      write (6,107) velr,veli
      write (8,101) Xvel,Yvel,velr,veli

      Go to 92

c-----------
      End If
c-----------

c-----
c Done
c-----

   99 Continue

      write (8,104) Null
      write (9,104) Null

      close (8)
      close (9)
      close (3)

 100  Format (1x,i3,20(1x,f15.10))
 101  Format (20(1x,f7.3))
 102  Format (1x,i3,20(1x,f8.5))
 104  Format (1x,i3,20(1x,f9.5))
 106  Format (" Real and imaginary parts of integral of the ",/
     +       ," boundary traction :",2(1x,f10.6))
 107  Format (" Real and imaginary parts of the velocity:",/
     +        ,2(1x,f10.5))
 150  format (1X,10(1X,f10.5))

      Stop
      End
