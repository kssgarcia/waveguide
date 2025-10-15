      program sgf_ax_1p_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c==============================================
c Driver for the periodic
c axisymmetric Green's function of Stokes flow
c
c Symbols
c ------
c
c  RL: Period
c
c  Iopt = 1  only the Green's function
c      ne 1  Green's function and stress tensor
c
c==============================================

      Implicit Double Precision (A-H,O-Z)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi
      pi8 = 8.0D0*pi

c--------
c prepare
c-------

      Iopt = 2      ! will compute all Green's functions

 98   Continue

      write (6,*) 
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 for one evaluation"
      write (6,*) "       2 to test the integral identities"
      write (6,*) "          on a circle "
      write (6,*) "       3 to test the integral identities"
      write (6,*) "         on a periodic line"
      write (6,*) "-----------------------------------------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

      If(menu.eq.2) then
       open (3,file="sgf_ax_1p.dat")
       read (3,*) RL
      Else
       write (6,*)
       write (6,*) " Please enter the period: L"
       write (6,*) "                 0 to quit "
       write (6,*) "---------------------------"
       read  (5,*) RL
      End If

      If(RL.eq.0) Go to 99

      If(menu.eq.2) then
       read (3,*) Nsum,Np
      Else
       write (6,*)
       write (6,*) " Please enter Nsum and Np "
       write (6,*) " for the summation"
       write (6,*) "            0 to quit "
       write (6,*) "----------------------"
       read (5,*) Nsum,Np
      End If

      If(Nsum.eq.0) go to 99
      If(Np  .eq.0) go to 99

c---------------
c one evaluation
c---------------

      If(menu.eq.1) then

      write (6,*)
      write (6,*) " Please enter x0 and s0 "
      write (6,*) "  99 for either to quit "
      write (6,*) "-------------------------"
      read (5,*) x0,s0

      If(x0.eq.99) Go to 99
      If(s0.eq.99) Go to 99

 97   Continue

      write (6,*)
      write (6,*) " Please enter x and s "
      write (6,*) " 99 for either to quit "
      write (6,*) "-----------------------"
      read (5,*) x,s

      If(x.eq.99) Go to 99
      If(s.eq.99) Go to 99

      call sgf_ax_1p
     +
     +   (Iopt
     +   ,x0,s0
     +   ,x,s
     +   ,RL
     +   ,Nsum,Np
     +   ,SXX,SXY
     +   ,SYX,SYY
     +   ,QXXX,QXXY,QXYX,QXYY
     +   ,QYXX,QYXY,QYYX,QYYY
     +   ,PXX,PXY,PYX,PYY
     +   ,Iaxis
     +   )

      write (6,*)
      write (6,*) " Velocity at (x,s) due to a line ring at (x0,s0)"
      write (6,*)
      write (6,*) " Axial point forces:"
      write (6,*) " -------------------"
      write (6,100) Sxx,Syx

      write (6,*)
      write (6,*) " Radial point forces:"
      write (6,*) " --------------------"
      write (6,100) Sxy,Syy

      If(Iopt.ne.1) then

       write (6,*)
       write (6,*) " Velocity at (x,s) due to a stresslet at (x0,s0)"
       write (6,*)
       write (6,*) " Stresslet strength xx:"
       write (6,100) Qxxx,Qyxx
       write (6,*) " Stresslet strength xs:"
       write (6,100) Qxxy,Qyxy
       write (6,*) " Stresslet strength sy:"
       write (6,100) Qxyy,Qyyy
       write (6,*) " Stresslet strength sy:"
       write (6,100) Qxyy,Qyyy
       write (6,*)

      End If
     
      Go to 97

c---------
      Else
c---------

c========================================
c Test of integral identities on a circle
c and a sinusoidal line
c========================================

 95   Continue

      write (6,*) " Testing identities (4), (5), (8), (9), (11)"
      write (6,*) 

      If(menu.eq.2) then
       read (3,*)
       read (3,*) x0,s0
      Else
       write (6,*) " Please enter x0 and s0"
       write (6,*) " 99 for one to quit "
       write (6,*) "--------------------"
       read (5,*) x0,s0
      End If

      If(x0.eq.99) Go to 99
      If(s0.eq.99) Go to 99

      If(menu.eq.2) then
       read  (3,*) mint
      Else
       write (6,*) " Enter the number of integration points"
       write (6,*) "  0 to quit"
       write (6,*) "-----------"
       read  (5,*) mint
      End If

      If(mint.eq.0) Go to 99

c----
      If(menu.eq.2) then
c----

c     write (6,*)
c     write (6,*) " Please enter the x and s coordinates"
c     write (6,*) "        of the center of the test circle   "
c     write (6,*) "        99 for either to quit "
c     write (6,*) "-------------------------------------------"
c     read (5,*) xcnt,scnt
      read (3,*) xcnt,scnt
      If(xcnt.eq.99) Go to 99
      If(scnt.eq.99) Go to 99

c     write (6,*)
c     write (6,*) " Please enter the radius of the test circle"
c     write (6,*) " 0 to quit"
c     write (6,*) "------------------------------"
c     read  (5,*) rad
      read  (3,*) rad
      If(rad.eq.0) Go to 99

c----
      Else If(menu.eq.3) then
c----

      write (6,*)
      write (6,*) " Enter the sigma position of the line"
      write (6,*) "    and amplitude of the sinusoid"
      write (6,*) " --------------------------------"
      read  (5,*) sline,amp

c----
      End If
c----

c---
c will verify that
c
c     SUM1 = 0
c     SUM2 = 0
c     SUM3 = 0 or 1
c     SUM4 = 0
c     SUM5 = 0
c     SUM6 = 0 or 1
c---

      If(menu.eq.2) then
        Dth = pi2/mint
      Else If(menu.eq.3) then
        wn = pi2/RL    ! wave number
        dx = RL/mint
      End If

      sum1 = 0.0D0
      sum2 = 0.0D0
      sum3 = 0.0D0
      sum4 = 0.0D0
      sum5 = 0.0D0
      sum6 = 0.0D0

      Do i=1,mint

        If(menu.eq.2) then
          th  = (i-0.5D0)*Dth
          vnx = Dcos(th)          ! unit normal
          vns = Dsin(th)
          x = xcnt + rad*vnx
          s = scnt + rad*vns
        Else If(menu.eq.3) then
          vnx = -amp*wn*cos(wn*x0) ! normal vector away from axis
          vns =  1.0D0
          x = (i-1.0D0)*dx+0.01D0
          s = sline + amp*sin(wn*x)
        End If

        call sgf_ax_1p
     +
     +    (Iopt
     +    ,x0,s0
     +    ,x,s
     +    ,RL
     +    ,Nsum,Np
     +    ,SXX,SXY
     +    ,SYX,SYY
     +    ,QXXX,QXXY,QXYX,QXYY
     +    ,QYXX,QYXY,QYYX,QYYY
     +    ,PXX,PXY,PYX,PYY
     +    ,Iaxis
     +    )

        sum1 = sum1 + (SXX *vnx + SYX *vns)*s
        sum2 = sum2 + (SXY *vnx + SYY *vns)*s

        sum3 = sum3 + (Qxxx*vnx + Qyxx*vns)*s
        sum4 = sum4 + (Qxxy*vnx + Qyxy*vns)*s
        sum5 = sum5 + (Qxyx*vnx + Qyyx*vns)*s
        sum6 = sum6 + (Qxyy*vnx + Qyyy*vns)*s

c       write (6,101) i,SUM1,SUM2

      End Do

      If(menu.eq.2) fc = rad*Dth/(pi8*s0)
      If(menu.eq.3) fc = dx/(pi8*s0)

      SUM1 = SUM1 * fc
      SUM2 = SUM2 * fc
      SUM3 = SUM3 * fc
      SUM4 = SUM4 * fc
      SUM5 = SUM5 * fc
      SUM6 = SUM6 * fc

      write (6,*) " -----------------------------------"
      write (6,*) " Identity (4): These should be: 0, 0:"
      write (6,100) sum1,sum2
      write (6,*) " -----------------------------------"
      write (6,*) " Identity (5): These should be:  0(1),  0"
      write (6,*) "                                  0   0(1)"
      write (6,100) sum3,sum4
      write (6,100) sum5,sum6

c----------------------------------
c Second set of integral identities
c----------------------------------

c---
c will verify that
c
c  SUM1 = 0
c  SUM2 = 0
c  SUM3 = 0 or -1
c  SUM4 = can be anything
c  SUM5 = 0 or -1
c  SUM6 = can be anything
c  SUM7 = 0
c  SUM8 = 0 or -1
c---

      SUM1 = 0.0D0
      SUM2 = 0.0D0
      SUM3 = 0.0D0
      SUM4 = 0.0D0
      SUM5 = 0.0D0
      SUM6 = 0.0D0
      SUM7 = 0.0D0
      SUM8 = 0.0D0

      Do i=1,mint

        If(menu.eq.2) then
          tht  = (i-0.5D0)*Dth
          vnx = Dcos(tht)          ! unit normal
          vns = Dsin(tht)
          x = xcnt + rad*vnx
          s = scnt + rad*vns
        Else If(menu.eq.3) then
          vnx = -amp*wn*cos(wn*x0) ! normal vector away from axis
          vns =  1.0D0
          x = (i-1.0D0)*dx+0.01D0
          s = sline + amp * dsin(wn*x)
        End If

        call sgf_ax_1p
     +
     +    (Iopt
     +    ,x,s
     +    ,x0,s0
     +    ,RL
     +    ,Nsum,Np
     +    ,SXX,SXY
     +    ,SYX,SYY
     +    ,QXXX,QXXY,QXYX,QXYY
     +    ,QYXX,QYXY,QYYX,QYYY
     +    ,PXX,PXY,PYX,PYY
     +    ,Iaxis
     +    )

        SUM1 = SUM1 + SXX *vnx + SXY *vns
        SUM2 = SUM2 + SYX *vnx + SYY *vns

        SUM3 = SUM3 + QXXX*vnx + QXXY*vns
        SUM4 = SUM4 + QXYX*vnx + QXYY*vns
        SUM5 = SUM5 + QYXX*vnx + QYXY*vns
        SUM6 = SUM6 + QYYX*vnx + QYYY*vns

        SUM7 = SUM7 + PXX *vnx + PXY *vns
        SUM8 = SUM8 + PYX *vnx + PYY *vns

      End Do

      If(menu.eq.2) fc = rad*Dth/pi8
      If(menu.eq.3) fc = dx/pi8

      SUM1 = SUM1 * fc
      SUM2 = SUM2 * fc
      SUM3 = SUM3 * fc
      SUM4 = SUM4 * fc
      SUM5 = SUM5 * fc
      SUM6 = SUM6 * fc
      SUM7 = SUM7 * fc
      SUM8 = SUM8 * fc

      write (6,*) " -----------------------------------"
      write (6,*) " Identity (8): These should be 0, 0:"
      write (6,100) SUM1,SUM2
      write (6,*) " ------------------------------------------"
      write (6,*) " Identity (9): These should be: 0(-1), anything:"
      write (6,100) sum3,sum4
      write (6,*) " ------------------------------------------"
      write (6,*) " Identity (9): These should be 0, anything:"
      write (6,100) sum5,sum6
      write (6,*) " ---------------------------------------"
      write (6,*) " Identity (11): These should be 0, 0(-1):"
      write (6,100) sum7,sum8
      write (6,*) " ---------------------------------------"

      If(menu.ne.2) Go to 95

c-----------
      End If
c-----------

      close (3)
 
      Go to 98    ! return to the main menu

c-----
c Done
c-----

 99   Continue

 100  Format (3(1x,f15.8))
 101  Format (1x,i3,8(1x,f8.5))

      Stop
      End
