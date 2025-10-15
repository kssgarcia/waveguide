      program sgf_ax_w_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c FDLIB BEMLIB
c
c Driver for the axisymmetric Green's function
c of Stokes flow in a semi-infinite domain bounded
c by a plane wall located at x = wall
c
c  Iopt = 1  only the Green's function
c      ne 1  Green's function and stress tensor
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi
      pi8 = 8.0D0*pi

c--------
c prepare
c--------

      Iopt = 2      ! will compute all Green's functions
      wall = 0.0;

 98   Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for one evaluation"
      write (6,*) " 2 to test the integral identities "
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99
      if(menu.gt.2) Go to 98

c---------------
c one evaluation
c---------------

      if(menu.eq.1) then

      write (6,*)
      write (6,*) " Please enter the x location of the wall"
      write (6,*) " 99 to quit "
      write (6,*) "------------"
      read (5,*) wall

      if(wall.eq.99) Go to 99

      write (6,*)
      write (6,*) " Please enter x0 and s0 "
      write (6,*) " 99 for either to quit "
      write (6,*) "-------------------------"
      read (5,*) x0,s0

      if(x0.eq.99) Go to 99
      if(s0.eq.99) Go to 99

 97   Continue

      write (6,*)
      write (6,*) " Please enter x and s "
      write (6,*) " 99 for either to quit "
      write (6,*) "-----------------------"
      read (5,*) x,s

      If(x.eq.99) Go to 99
      If(s.eq.99) Go to 99

      call sgf_ax_w
     +
     +  (Iopt
     +  ,x0,s0
     +  ,x,s
     +  ,wall
     +  ,Sxx,Sxy
     +  ,Syx,Syy
     +  ,QXXX,QXXY,QXYX,QXYY
     +  ,QYXX,QYXY,QYYX,QYYY
     +  ,PXX,PXY,PYX,PYY
     +  ,Iaxis
     +  )

      write (6,*)
      write (6,*) " Axial point forces:"
      write (6,*) " -------------------"
      write (6,100) Sxx,Syx

      write (6,*)
      write (6,*) " Radial point forces:"
      write (6,*) " --------------------"
      write (6,100) Sxy,Syy

      if(Iopt.ne.1) then

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

      end if
     
      Go to 97

      end if

c----------------------------
      if(menu.eq.2) then
c----------------------------

      write (6,*) " Testing integral identities"

c----------------------------------------
c test of integral identities on a circle
c----------------------------------------

 95   Continue

c     open (3,file="sgf_ax_w.dat")
c      read (3,*) wall
c      read (3,*) x0,s0
c      read (3,*) mint
c      read (3,*) xcnt,scnt
c      read (3,*) rad
c      close (3)


c     write (6,*)
c     write (6,*) " Enter x0 and s0 "
c     write (6,*)
c     write (6,*) " 99 for either to quit "
c     write (6,*) " ----------------------"
c     read (5,*) x0,s0
      x0 = 0.1
      s0 = 0.9

      if(x0.eq.99) Go to 99
      if(s0.eq.99) Go to 99

c     write (6,*)
c     write (6,*) " Enter the x and s coordinates"
c     write (6,*) " of the center of the test circle   "
c     write (6,*)
c     write (6,*) " 99 for either to quit "
c     write (6,*) " ----------------------"
c     read (5,*) xcnt,scnt
      xcnt = 0.21
      scnt = 1.18

      if(xcnt.eq.99) Go to 99
      if(scnt.eq.99) Go to 99

c     write (6,*)
c     write (6,*) " Enter the radius of the test circle"
c     write (6,*)
c     write (6,*) " 0 to quit"
c     write (6,*) "----------"
c     read  (5,*) rad
      rad = 0.38

      if(rad.eq.0) Go to 99

c     write (6,*) " Enter the number of integration points"
c     write (6,*) 
c     write (6,*) " 0 to quit"
c     write (6,*) "----------"
c     read  (5,*) mint
      mint = 256

      if(mint.eq.0) Go to 99

c---
c first set of identities
c---

      Dth = pi2/mint

      sum1 = 0.0D0
      sum2 = 0.0D0
      sum3 = 0.0D0
      sum4 = 0.0D0
      sum5 = 0.0D0
      sum6 = 0.0D0

      Do i=1,mint

        th  = (i-0.5D0)*Dth
        vnx = Dcos(th)          ! unit normal
        vns = Dsin(th)

        x = xcnt + rad*vnx
        s = scnt + rad*vns

        call sgf_ax_w
     +
     +    (Iopt
     +    ,x,s
     +    ,x0,s0
     +    ,wall
     +    ,Sxx,SXY
     +    ,SYX,SYY
     +    ,QXXX,QXXY,Qxyx,Qxyy
     +    ,QYXX,QYXY,Qyyx,Qyyy
     +    ,PXX,PXY,PYX,PYY
     +    ,Iaxis
     +    )

        sum1 = sum1 + (vnx*Sxx + vns*SYX)*s
        sum2 = sum2 + (vnx*SXY + vns*SYY)*s

        sum3 = sum3 + (vnx*Qxxx + vns*Qyxx)*s
        sum4 = sum4 + (vnx*Qxxy + vns*Qyxy)*s
        sum5 = sum5 + (vnx*Qxyx + vns*Qyyx)*s
        sum6 = sum6 + (vnx*Qxyy + vns*Qyyy)*s

c       write (6,101) i,SUM1,SUM2

      End Do

      fc = rad*Dth/(pi8*s0)

      SUM1 = SUM1 * fc
      SUM2 = SUM2 * fc
      SUM3 = SUM3 * fc
      SUM4 = SUM4 * fc
      SUM5 = SUM5 * fc
      SUM6 = SUM6 * fc

      write (6,*) " -----------------------------------"
      write (6,*) " These should be: 0, 0:"
      write (6,100) sum1,sum2
      write (6,*) " -----------------------------------"
      write (6,*) " These should be:  0(1),  0"
      write (6,*) "                   0   0(1)"
      write (6,100) sum3,sum4
      write (6,100) sum5,sum6
      write (6,*) " -----------------------------------"

c----------------------------------
c second set of integral identities
c----------------------------------

      SUM1 = 0.0D0
      SUM2 = 0.0D0
      SUM3 = 0.0D0
      SUM4 = 0.0D0
      SUM5 = 0.0D0
      SUM6 = 0.0D0
      sum7 = 0.0D0
      sum8 = 0.0D0

      Do i=1,mint

        th = (i-0.5D0)*Dth
        vnx = Dcos(th)
        vns = Dsin(th)
        x   = xcnt + rad*vnx
        s   = scnt + rad*vns

        call sgf_ax_w
     +
     +    (Iopt
     +    ,x0,s0
     +    ,x,s
     +    ,wall
     +    ,Sxx,SXY,SYX,SYY
     +    ,QXXX,QXXY,QXYX,QXYY
     +    ,QYXX,QYXY,QYYX,QYYY
     +    ,PXX,PXY,PYX,PYY
     +    ,Iaxis
     +    )

        SUM1 = SUM1 + Sxx *vnx + SXY *vns
        SUM2 = SUM2 + SYX *vnx + SYY *vns

        SUM3 = SUM3 + QXXX*vnx + QXXY*vns
        SUM4 = SUM4 + QXYX*vnx + QXYY*vns
        SUM5 = SUM5 + QYXX*vnx + QYXY*vns
        SUM6 = SUM6 + QYYX*vnx + QYYY*vns

        sum7 = sum7 + PXX *vnx + PXY *vns
        sum8 = sum8 + PYX *vnx + PYY *vns

      End Do

      cf = rad*dth/pi8

      SUM3 = SUM3 * cf
      SUM4 = SUM4 * cf
      SUM5 = SUM5 * cf
      SUM6 = SUM6 * cf
      sum7 = sum7 * cf
      sum8 = sum8 * cf

      write (6,*) " --------------------------"
      write (6,*) " Should be 0, 0:"
      write (6,100) sum1,sum2
      write (6,*) " ---------------------------"
      write (6,*) " Should be: 0(-1), anything:"
      write (6,100) sum3,sum4
      write (6,*) " ---------------------------"
      write (6,*) " Should be 0, anything:"
      write (6,100) sum5,sum6
      write (6,*) " ---------------------------"
      write (6,*) " Should be 0, 0(-1):"
      write (6,100) sum7,sum8
      write (6,*) " ---------------------------"

c     Go to 95

c-----------
      end if
c-----------

      Go to 98    ! return to the main menu

c-----
c done
c-----

 99   Continue

 100  format (3(1x,f15.8))
 101  format (1x,i3,8(1x,f8.5))

      stop
      end
