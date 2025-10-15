      program sgf_ax_1p_ct_dr

c--------------------------------------------
c FDLIB and BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c--------------------------------------------

c---------------------------------------------
c Driver program for the
c axisymmetric periodic Green's function
c of Stokes flow inside a circular cylinder
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension ux(500),flowx(500)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi
      pi8 = 8.0D0*pi

      Null = 0

c--------
c prepare
c--------

      open (2,file="data.out")

c-------
c launch
c-------

 98   Continue

      write (6,*)
      write (6,*) " Enter: "
      write (6,*)
      write (6,*) " 1 for one evaluation"
      write (6,*) " 2 to test the continuity equation"
      write (6,*) " 3 for profile - flow rate"
      write (6,*) " 4 to compute the drag on the wall"
      write (6,*) "   and the pressure drop over a period"
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read (5,*) menu

      if(menu.eq.0) Go to 99

c----------------------
c read input parameters
c----------------------

      open (1,file="sgf_ax_1p_ct.dat")

       read (1,*) sc
       read (1,*) RL
       read (1,*) Nsum,Np

      close (1)

c---------
      Go to 331    ! skip the interactive session
c---------

      write (6,*)
      write (6,*) "Enter the tube radius" 
      write (6,*) "---------------------"
      read  (5,*) sc

      write (6,*)
      write (6,*) "Enter the period" 
      write (6,*) "----------------"
      read  (5,*) RL

      write (6,*)
      write (6,*) "Enter Nsum and Np"
      write (6,*)
      write (6,*) "0 to quit"
      write (6,*) "----------"
c     read  (5,*) Nsum,Np
      Nsum = 5
      Np = 2

      If(Nsum.eq.0.or.Np.eq.1) Go to 99

 331  Continue

 93   Continue

      write (6,*) "Enter x0 and s0"
      write (6,*) "99 for either to quit "
      write (6,*) "----------------------"
      read (5,*) x0,s0

      If(x0.eq.99) Go to 99
      If(s0.eq.99) Go to 99

c---------------
c one evaluation
c---------------

      If(menu.eq.1) then

 97   Continue

      write (6,*) "Enter x and s"
      write (6,*) "99 for either to quit "
      write (6,*) "----------------------"
      read (5,*) x,s

      If(x.eq.99) Go to 99
      If(s.eq.99) Go to 99

      call sgf_ax_1p_ct 
     +
     +    (x,s 
     +    ,x0,s0   ! singularity
     +    ,sc
     +    ,RL
     +    ,Nsum,Np
     +    ,Gxx,Gxs
     +    ,Gsx,Gss
     +    )

      write (6,*)
      write (6,*) " Velocity at (x,s) due to a line ring at (x0,s0)"
      write (6,*)
      write (6,*) " Axial point forces:"
      write (6,*) " -------------------"
      write (6,100) Gxx,Gsx

      write (6,*)
      write (6,*) " Radial point forces:"
      write (6,*) " --------------------"
      write (6,100) Gxs,Gss

c---
c also compute the free-space GF
c for comparison
c----

      Go to 333

      Iopt = 1

      call sgf_ax_fs
     +
     +     (Iopt
     +     ,x,s
     +     ,x0,s0
     +     ,SXX,SXY
     +     ,SYX,SYY
     +     ,QXXX,QXXY,QXYX,QXYY
     +     ,QYXX,QYXY,QYYX,QYYY
     +     ,PXX,PXY,PYX,PYY
     +     ,Iaxis
     +     )

      write (6,*)
      write (6,*) " Velocity at (x,s) due to a line ring at (x0,s0)"
      write (6,*)
      write (6,*) " Axial point forces:"
      write (6,*) " -------------------"
      write (6,100) Sxx,Sxy

      write (6,*)
      write (6,*) " Radial point forces:"
      write (6,*) " --------------------"
      write (6,100) Syx,Syy

 333  Continue

      Go to 97

c----------------------------
      elseif(menu.eq.2) then
c----------------------------

 95   Continue

      write (6,*) "Enter the x and s coordinates"
      write (6,*) "of the center of the test circle   "
      write (6,*) "99 for either to quit "
      write (6,*) "----------------------"
      read (5,*) xcnt,scnt

      If(xcnt.eq.99) Go to 99
      If(scnt.eq.99) Go to 99

      write (6,*)
      write (6,*) "Please enter the radius of the test circle"
      write (6,*) "0 to quit"
      write (6,*) "---------"
      read  (5,*) rad

      If(rad.eq.0) Go to 99

      write (6,*) "Please enter the number of integration points"
      write (6,*) "0 to quit"
      write (6,*) "---------"
      read  (5,*) mint

      If(mint.eq.0) Go to 99

c---
c will verify that
c
c   SUM1 = 0
c   SUM2 = 0
c---

      Dth = pi2/mint

      SUM1 = 0.0D0
      SUM2 = 0.0D0

      Do i=1,mint

        th = (i-0.50D0)*Dth
        vnx = Dcos(th)       ! unit normal vector
        vns = Dsin(th)

        x = xcnt + rad*vnx
        s = scnt + rad*vns

        call sgf_ax_1p_ct 
     +
     +     (x,s
     +     ,x0,s0     ! singularity
     +     ,sc
     +     ,RL
     +     ,Nsum,Np
     +     ,Gxx,Gxy
     +     ,Gyx,Gyy
     +     )

        SUM1 = SUM1 + (GXX*vnx + GYX*vns)*s
        SUM2 = SUM2 + (GXY*vnx + GYY*vns)*s

      End Do

      write (6,*) " -----------------------------------"
      write (6,*) " Identity (4): These should be: 0, 0:"
      write (6,100) sum1,sum2


c---
c will verify that
c
c             SUM1 = 0
c             SUM2 = 0
c---

      write (6,*)
      write (6,*) " Second test of the continuity equation"
      write (6,*) " -------------------------------------"

      SUM1 = 0.0D0
      SUM2 = 0.0D0

      Do i=1,mint

        th  = (i-0.5D0)*Dth
        vnx = Dcos(th)          ! normal vector
        vns = Dsin(th)
        x   = xcnt + rad*vnx
        s   = scnt + rad*vns

        call sgf_ax_1p_ct
     +
     +      (x0,s0
     +      ,x,s
     +      ,sc
     +      ,RL
     +      ,Nsum,Np
     +      ,Sxx,Sxy
     +      ,Syx,Syy
     +      )

        SUM1 = SUM1 + SXX *vnx + SXY *vns
        SUM2 = SUM2 + SYX *vnx + SYY *vns

      End Do

      write (6,*) " -----------------------------------"
      write (6,*) " Identity (7): These should be 0, 0:"
      write (6,100) sum1,sum2

      Go to 95

c-----------------------------
      elseif (menu.eq.3) then
c-----------------------------

 96   Continue

      write (6,*)
      write (6,*) "Enter the x position of the profile" 
      write (6,*) "99 to quit"
      write (6,*) "----------"
      read  (5,*) x

      If(x.eq.99) Go to 99

      write (6,*)
      write (6,*) " How many profile points?"
      write (6,*)
      write (6,*) " Please choose an even number"
      write (6,*) " Enter 0 to quit"
      write (6,*) "----------------"
      read  (5,*) Nprf

      If(Nprf.eq.0) Go to 99

      Nprf1 = Nprf+1

      ds = sc/Nprf

      write (6,*) Nprf1
      write (2,*) Nprf1

c---
c generate the profile
c---

      Do i=1,Nprf1

        s = (i-1.0D0)*ds
        if(abs(s).le.0.001) s = 0.00001D0

        call sgf_ax_1p_ct 
     +
     +     (x,s
     +     ,x0,s0  ! singularity point
     +     ,sc
     +     ,RL
     +     ,Nsum,Np
     +     ,gxx,gxs
     +     ,gsx,gss
     +     )

           ux(i) = gxx
        flowx(i) = gxx*s

        write (6,101) i,s,ux(i),flowx(i)
        write (2,101) i,s,ux(i),flowx(i)

      End Do

c---
c Compute the flow rate by Simpson's rule
c---

      flowx(Nprf1) = 0.0D0

      Flow = 0.0D0

      Do i=2,Np,2
        Flow = Flow + flowx(i-1)+4.0D0*flowx(i)+flowx(i+1)
      End Do

      Flow = pi2*flow*ds/3.0D0

      write (6,102) Flow

      Go to 96

c-----------------------------
      elseif (menu.eq.4) then
c-----------------------------

 92   Continue

      write (6,*)
      write (6,*) "How many divisions in a period?"
      write (6,*)
      write (6,*) "Enter 0 to quit"
      write (6,*) "----------------"
      read  (5,*) Ndiv

      if(Ndiv.eq.0) Go to 99

      write (6,*)
      write (6,*) " Please enter epsilon for numerical differentiation"
      write (6,*) " of the shear stress"
      write (6,*) "--------------------"
      read  (5,*) eps

c-----
c viscosity
c-----

      visc = 1.345D0     ! this value is irrelevant

c------
c will extrapolate to eps = 0
c------

      Ipass = 1

  88  Continue

      If(Ipass.eq.2) eps = 0.5D0*eps

      s = sc-eps

      step = RL/Ndiv

      Drag = 0.0D0

      Do i=1,Ndiv

        x = 0.01234+(i-1.0D0)*step

        call sgf_ax_1p_ct 
     +
     +     (x,s
     +     ,x0,s0  ! singularity point
     +     ,sc
     +     ,RL
     +     ,Nsum,Np
     +     ,gxx,gxs
     +     ,gsx,gss
     +     )

        shear = (0.0D0-gxx)/eps

        Drag = Drag + shear

        write (6,101) i,x,s,shear

      End Do
      
      Drag = visc * (pi2*sc) * (Drag*step)  ! multiply by the wall area

      Drag = Drag/(pi8*visc)   ! definition of the point force involves
                               ! a factor of 8*pi

      If(Ipass.eq.1) then
        Drag1 = Drag
        Ipass = 2
        Go to 88
      Else
        Drag2 = Drag
      End If

      Drag = 2.0D0*Drag2-Drag1

c------------------------------------------- 
c compute the force due to the pressure drop
c-------------------------------------------

      Dp = - 4.0D0 * s0 * (sc*sc-s0*s0)/sc**4

      PD =  pi*sc*sc*Dp  ! multiply by the cross-sectional area

c---
c x-force balance
c---

      sum = PD+Drag

      sum = sum/(pi2*s0)  ! normalize

c---
c display
c---

      write (6,*) 
      write (6,745) Drag,PD,sum

      Go to 93
      Go to 92

c------------
      end if
c------------

  99  Continue

c-----
c Done
c-----

      write (2,100) Null
      close (2)

  100 Format (5(1x,f15.10))
  101 Format (1x,i3,5(1x,f15.10))
  102 Format ("Flow Rate = ",f15.10)
  103 Format (5(1x,f10.5))
  745 Format (    "Drag =", f15.10
     +        ,//,"PD   =", f15.10
     +        ,//,"sum  =", f15.10)

      Stop
      End
