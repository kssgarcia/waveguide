      program lgf_2d_ww_dr

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c==========================================

c---------------------------------------------------------
c Driver for the 2-D Green's function of Laplace's equation
c in a domain confined between two parallel plates
c
c Set Iopt  = 1 to compute G
c     Iopt ne 2 to compute G and gradient
c---------------------------------------------------------

      Implicit double precision (a-h,o-z)

      Dimension u1(900)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c--------
c prepare
c--------

      Iopt = 2   ! will need G and grad(G)

c---
c launch
c---

      write (6,*) 
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 for one value "
      write (6,*) " 2 to test integral identities "
      write (6,*) " 3 y-profile and flow rate"
      write (6,*) " --------------------------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

  93  Continue

      write (6,*)
      write (6,*) ' The first wall is located at y = wall1'
      write (6,*)
      write (6,*) ' Please enter: wall1'
      write (6,*) '          99 to quit'
      write (6,*) ' -------------------'
      read  (5,*) wall1

      write (6,*)
      write (6,*) ' The second wall is located at y = wall2'
      write (6,*)
      write (6,*) ' Please enter: wall2'
      write (6,*)
      write (6,*) ' 99 to quit'
      write (6,*) ' ----------'
      read  (5,*) wall2

      If(wall2.lt.wall1) then
       write (6,*) 
       write (6,*) " Walls are poorly defined"
       write (6,*) " Must be wall2 > wall1"
       Go to 93
      End If

      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 0 to quit'
      write (6,*) ' 1 for the Green function'
      write (6,*) ' 2 for the Neumann function'
      write (6,*) '---------------------------'
      read  (5,*) Ign

 98   Continue

      write (6,*)
      write (6,*) " Please enter the singularity point: x0, y0"
      write (6,*)
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
      read  (5,*) x0,y0

      If(x0.eq.99) Go to 99
      If(y0.eq.99) Go to 99

c-----------------------
      If(menu.eq.1) then
c-----------------------

c---------------
c One evaluation
c---------------


      write (6,*)
      write (6,*) " Enter the field-point coordinates: x, y"
      write (6,*)
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
      read  (5,*) x,y

      If(x.eq.99) Go to 99
      If(y.eq.99) Go to 99

      call lgf_2d_ww 
     +
     +  (Iopt
     +  ,Ign
     +  ,x,y
     +  ,x0,y0
     +  ,wall1,wall2
     +  ,G
     +  ,Gx,Gy
     +  )

      write (6,*) ' ---------------------------'
      write (6,*) ' Green function and gradient'
      write (6,*) ' ---------------------------'

      write (6,100) G
      write (6,100) Gx,Gy
      write (6,100) 

c----------------------------
      Else If(menu.eq.2) then
c----------------------------

c-----------------------------
c Test of integral identities
c
c sm1 should be 0 or -1
c-----------------------------


      write (6,*)
      write (6,*) " Please enter:"
      write (6,*)
      write (6,*) "  the center of the test circle "
      write (6,*) "  the radius of  the test circle"
      write (6,*) "  the number of integration points "
      write (6,*) " ----------------------------------"
      read  (5,*) xcnt,ycnt,rad,mint

      mint1 = ming+1

      dth = pi2/(mint1-1.0D0)

      sm1 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = Dcos(th)
        sn = Dsin(th)
        x  = xcnt + rad*cs
        y  = ycnt + rad*sn

        call lgf_2d_ww
     +
     +   (Iopt
     +   ,Ign
     +   ,x,y
     +   ,x0,y0
     +   ,wall1,wall2
     +   ,G
     +   ,Gx,Gy
     +   )

        vnx = cs  ! normal vector pointing outward
        vny = sn

        sm1 = sm1 + Gx*vnx + Gy*vny

      End Do

      cf  = dth*rad
      sm1 = sm1 * cf

      write (6,100)
      write (6,100) sm1
      write (6,100)

c-----------------------------
      Else If (menu.eq.3) then
c-----------------------------

c------------------------
c y-profile of x-velocity
c and flow rate
c------------------------

      write (6,*)
      write (6,*) " Please enter: "
      write (6,*)
      write (6,*) " the x position"
      write (6,*) " the number of profile points (should be even)"
      write (6,*) " ---------------------------------------------"
      read  (5,*) x,mprof

      mprof1 = mprof+1

      dy = (wall2-wall1)/(mprof1-1.0D0)

      Do i=1,mprof1

        y = wall1+(i-1.0D0)*dy

        call lgf_2d_ww 
     +
     +    (Iopt
     +    ,Ign
     +    ,x,y
     +    ,x0,y0
     +    ,wall1,wall2
     +    ,G
     +    ,Gx,Gy
     +    )

        u1(i) = Gx

        write (6,101) i,y,u1(i)

      End Do

      Flow = 0.0D0

      Do i=1,mprof
        Flow = Flow + u1(i) + u1(i+1)
      End Do

      Flow = 0.5D0* dy * Flow

      write (6,102) Flow

c-----------
      End If
c-----------

      Go to 98

  99  Continue

c-----
C Done
c-----

 100  Format (10(1x,f12.6))
 101  Format (1x,i3,4(1x,f15.10))
 102  Format (" Flow rate =",f20.10)

      Stop
      End
