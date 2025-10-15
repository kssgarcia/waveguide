      program lgf_2d_crc_dr

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
c Driver for the 2-D Neumann function of Laplace's equation
c in the exterior of a circle
c
c Set Iopt  = 1 to compute G
c     Iopt ne 2 to compute G and gradient
c---------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pi2 = 2.0D0 * pi
      pi4 = 4.0D0 * pi

c--------
c prepare
c--------

      Iopt = 2

c---
c launch
c---

 98   continue

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 for one value "
      write (6,*) " 2 to test integral identities "
      write (6,*) " ------------------------------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Please enter:"
      write (6,*) 
      write (6,*) " The coordinates of the circle center"
      write (6,*) "                    the circle radius"
      write (6,*) " -------------------------------------"
      read  (5,*) xc,yc
      read  (5,*) a

      write (6,*) 
      write (6,*) " Please enter the coordinates"
      write (6,*) " of the singular point: x0,y0"
      write (6,*) 
      write (6,*) " 99 for either to quit"
      write (6,*) " ---------------------"
      read  (5,*) x0,y0

      If(x0.eq.99) Go to 99
      If(y0.eq.99) Go to 99

c---------------
c One evaluation
c---------------

      If(menu.eq.1) then

      write (6,*) 
      write (6,*) " Enter field-point coordinates x, y"
      write (6,*) 
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read  (5,*) x,y

      If(x.eq.99) Go to 99
      If(y.eq.99) Go to 99

      call lgf_2d_crc
     +
     +  (Iopt
     +  ,x,y
     +  ,x0,y0
     +  ,xc,yc
     +  ,a
     +  ,G
     +  ,Gx,Gy
     +  )


      write (6,*)  " ---------------------------"
      write (6,*)  " Green function and gradient"
      write (6,100) G
      write (6,100) Gx,Gy
      write (6,*)  " ---------------------------"

c-----------------------------
c  Test of integral identities
c-----------------------------

      Else If(menu.eq.2) then


      write (6,*)
      write (6,*) " Enter the center of the test circle "
      write (6,*) " -----------------------------------"
      read (5,*) xcnt,ycnt

      write (6,*)
      write (6,*) " Enter the radius of the test circle"
      write (6,*) " -----------------------------------"
      read (5,*) rad

      write (6,*)
      write (6,*) " Enter the number of integration points "
      write (6,*) " --------------------------------------"
      read (5,*) mint

      mint1  = mint+1

      dth = pi2/(mint1-1.0D0)

      sm1 = 0.0D0

      Do i=1,mint

        th = (i-1.0D0)*dth
        cs = Dcos(th)
        sn = Dsin(th)
        x  = xcnt + rad * cs
        y  = ycnt + rad * sn

        call lgf_2d_crc
     +
     +    (Iopt
     +    ,x,y
     +    ,x0,y0
     +    ,xc,yc
     +    ,a
     +    ,G
     +    ,Gx,Gy
     +    )

        vnx = cs   ! normal vector
        vny = sn

        sm1 = sm1 + Gx*vnx + Gy*vny

        write (6,*) sm1

      End Do

        write (6,*) sm1
      cf = dth*rad

      sm1 = sm1 * cf

      write (6,*)
      write (6,*) " Should be 0 or -1.0"
      write (6,*)
      write (6,100) sm1
      write (6,100)

c---
      End If
c---

      Go to 98

c-----
c Done
c-----

  99  Continue

 100  Format (10(1x,f12.6))
 101  Format (1x,i3,4(1x,f15.10))
 102  Format (" Flow rate =",f20.10)

      Stop
      End
