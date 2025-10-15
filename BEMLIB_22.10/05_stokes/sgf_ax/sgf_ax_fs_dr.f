      program sgf_ax_fs_dr

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
c Driver for the axisymmetric Green's function
c of Stokes flow in free space
c 
c Iopt = 1  compute only the Green's function
c     ne 1  compute the Green's function
c           and the stress tensor
c--------------------------------------------

      Implicit Double Precision (A-H,O-Z)

c----------
c constants
c----------

      pi  = 3.1415 92653 58979 32384 D0
      pi2 = 2.0D0*pi
      pi8 = 8.0D0*pi

c--------
c prepare
c--------

      Iopt = 2      ! will compute G, p, Q

 98   Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for one evaluation"
      write (6,*) " 2 to test integral identities"
      write (6,*) " 0 to quit"
      write (6,*) "----------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c---------------
c one evaluation
c---------------

      if(menu.eq.1) then

c     write (6,*)
c     write (6,*) " Enter x0 and s0"  
c     write (6,*)
c     write (6,*) " 99 for either to quit"
c     write (6,*) " ---------------------"
c     read (5,*) x0,s0
      x0 = 0.0
      s0 = 0.9

      if(x0.eq.99) Go to 99
      if(s0.eq.99) Go to 99

 97   Continue

c     write (6,*)
c     write (6,*) " Enter x and s"
c     write (6,*)
c     write (6,*) " 99 for either to quit"
c     write (6,*) " ----------------------"
c     read (5,*) x,s
      x = 0.34
      s = 1.48

      if(x.eq.99) Go to 99
      If(s.eq.99) Go to 99

      call sgf_ax_fs
     +
     +   (Iopt
     +   ,x0,s0 
     +   ,x,s      ! singular point
     +   ,Sxx,Sxy
     +   ,Syx,Syy
     +   ,QXXX,QXXY,QXYX,QXYY
     +   ,QYXX,QYXY,QYYX,QYYY
     +   ,PXX,PXY,PYX,PYY
     +   ,Iaxis
     +   )

      write (6,*)
      write (6,*) " Velocity at (x0,s0) due to a line ring at (x,s)"
      write (6,*)
      write (6,*) " axial point forces:"
      write (6,*) " -------------------"
      write (6,100) Sxx,Syx

      write (6,*)
      write (6,*) " Radial point forces:"
      write (6,*) " --------------------"
      write (6,100) Sxy,Syy

      if(Iopt.ne.1) then

       write (6,*)
       write (6,*) " Velocity at (x0,s0) due to a stresslet at (x,s)"
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

c----------------------------
      else if(menu.eq.2) then
c----------------------------

c----------------------------------------
c test of integral identities on a circle
c----------------------------------------

 95   Continue

c     write (6,*) " Enter x0 and s0"
c     write (6,*)
c     write (6,*) " 99 for either to quit"
c     write (6,*) " ---------------------"
c     read (5,*) x0,s0
      x0 = 0.0
      s0 = 0.9

      if(x0.eq.99) Go to 99
      if(s0.eq.99) Go to 99

c     write (6,*)
c     write (6,*) " Enter the x and s coordinates"
c     write (6,*) " of the center of the test circle "
c     write (6,*)
c     write (6,*) " 99 for either to quit "
c     write (6,*) "-----------------------"
c     read (5,*) xcnt,scnt
      xcnt = 0.21
      scnt = 1.48

      if(xcnt.eq.99) Go to 99
      if(scnt.eq.99) Go to 99

c     write (6,*)
c     write (6,*) " Enter the radius of the test circle"
c     write (6,*)
c     write (6,*) " 0 to quit"
c     write (6,*) "----------"
c     read  (5,*) rad
      rad = 0.84

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

        call sgf_ax_fs
     +
     +    (Iopt
     +    ,x,s
     +    ,x0,s0
     +    ,SXX,SXY
     +    ,SYX,SYY
     +    ,QXXX,QXXY,QXYX,QXYY
     +    ,QYXX,QYXY,QYYX,QYYY
     +    ,PXX,PXY,PYX,PYY
     +    ,Iaxis
     +    )

        sum1 = sum1 + (vnx*SXX + vns*SYX)*s
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

        th = (i-0.5D0)*Dth
        vnx = Dcos(th)
        vns = Dsin(th)
        x   = xcnt + rad*vnx
        s   = scnt + rad*vns

        call sgf_ax_fs 
     +
     +     (Iopt
     +     ,x0,s0
     +     ,x,s
     +     ,SXX,SXY,SYX,SYY
     +     ,QXXX,QXXY,QXYX,QXYY
     +     ,QYXX,QYXY,QYYX,QYYY
     +     ,PXX,PXY,PYX,PYY
     +     ,Iaxis
     +     )

        SUM1 = SUM1 + SXX*vnx + SXY*vns
        SUM2 = SUM2 + SYX*vnx + SYY*vns

        SUM3 = SUM3 + QXXX*vnx + QXXY*vns
        SUM4 = SUM4 + QXYX*vnx + QXYY*vns

        SUM5 = SUM5 + QYXX*vnx + QYXY*vns
        SUM6 = SUM6 + QYYX*vnx + QYYY*vns

        SUM7 = SUM7 +  PXX*vnx +  PXY*vns
        SUM8 = SUM8 +  PYX*vnx +  PYY*vns

      End Do

      cf = rad*dth/pi8

      SUM3 = SUM3 * cf
      SUM4 = SUM4 * cf
      SUM5 = SUM5 * cf
      SUM6 = SUM6 * cf
      SUM7 = SUM7 * cf
      SUM8 = SUM8 * cf

      write (6,*) " ---------------------------"
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
      End If
c-----------

      Go to 98    ! return to the main menu

c-----
c done
c-----

 99   Continue

 100  Format (3(1x,f15.8))
 101  Format (1x,i3,8(1x,f8.5))

      Stop
      End
