      program lgf_2d_www_dr

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c---------------------------------------------------------
c Driver for the 2-D Neumann function of Laplace's equation
c in a domain confined between two parallel plates
c and another plate intersecting them at a right angle
c
c The first  plate is located at x = wall1
c The second plate is located at x = wall2
c The third  plate is located at y = wall3
c
c Set Iopt  = 1 to compute G
c     Iopt ne 2 to compute G and gradient
c---------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

c--------
c prepare
c--------

      Iopt = 2

c---
c launch
c---

 98   Continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 1 for one value "
      write (6,*) " 2 to test the integral identity "
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) menu

      if(menu.eq.0) Go to 99

 78   Continue

      write (6,*)
      write (6,*) " Please enter: wall1, wall2, wall3"
      write (6,*) " ---------------------------------"
      read  (5,*) wall1,wall2,wall3

      if(wall2.lt.wall1) then
       write (6,*) " Walls are defined improperly"
       write (6,*) " Must be wall2 >= wall1"
       Go to 78
      end if

      write (6,*)
      write (6,*) " Please enter the singularity point: x0,y0"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) x0,y0

      if(x0.eq.99) Go to 99
      if(y0.eq.99) Go to 99

c---------------
c One evaluation
c---------------

      if(menu.eq.1) then

      write (6,*)
      write (6,*) " Enter the field-point coordinates x, y"
      write (6,*)
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
      read  (5,*) x,y

      if(x.eq.99) Go to 99
      if(y.eq.99) Go to 99

      call lgf_2d_www
     +
     +   (Iopt
     +   ,x,y
     +   ,x0,y0
     +   ,wall1,wall2,wall3
     +   ,G
     +   ,Gx,Gy
     +   )

      write (6,100) 
      write (6,100) G
      write (6,100) Gx,Gy
      write (6,100) 
      write (6,100)

c-------------------------
c test integral identities
c-------------------------

      else if(menu.eq.2) then

      write (6,*)
      write (6,*) " Enter: "
      write (6,*)
      write (6,*) " the (x, y) coord of the center of the test circle"
      write (6,*) " the radius of  the test circle"
      write (6,*) " the number of integration points"
      write (6,*) " ---------------------------------"
      read (5,*) xcnt,ycnt,rad,mint

      mint1  = mint+1
      dth = pi2/(mint1-1.0D0)

      sm1 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = Dcos(th)
        sn = Dsin(th)
        x  = xcnt + rad*cs
        y  = ycnt + rad*sn

        call lgf_2d_www
     +
     +    (Iopt
     +    ,x,y
     +    ,x0,y0
     +    ,wall1,wall2,wall3
     +    ,G
     +    ,Gx,Gy
     +    )

        vnx = cs   ! normal vector pointing outward
        vny = sn

        sm1 = sm1 + Gx*vnx + Gy*vny

      End Do

      cf = dth*rad

      sm1 = sm1 * cf

      write (6,*) " should be -1:"
      write (6,100) sm1
      write (6,*)

c---
      end if
c---

      Go to 98 ! to repeat

c-----
c Done
c-----

  99  Continue

 100  Format (10(1x,f12.6))
 101  Format (1x,i3,4(1x,f15.10))
 102  Format (" Flow rate =",f20.10)

      Stop
      End
