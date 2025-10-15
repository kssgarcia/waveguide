      program lgf_2d_fs_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------------
c Driver for Green's function of Laplace's equation
c in free space
c--------------------------------------------------
 
      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159 265358D0
      pi2 = 2.0*pi
      pi4 = 4.0*pi

c-----------
c initialize
c-----------

      Iopt = 2   ! need G and gradient

c-------
c launch
c-------

  98  Continue

      write (6,*)
      write (6,*) '         MENU'
      write (6,*)
      write (6,*) ' Enter 0 to quit'
      write (6,*) '       1 for one evaluation'
      write (6,*) '       2 to test an integral identity'
      write (6,*) '-------------------------------------'
      read  (5,*) menu

      If(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) ' Enter the coordinates of the singular point'
      write (6,*) '                                  99 to quit'
      write (6,*) ' -------------------------------------------'
      read  (5,*) x0,y0

      If(x0.eq.99) Go to 99
      If(y0.eq.99) Go to 99

c---------------
c One evaluation
c---------------

      If(menu.eq.1) then

      write (6,*)
      write (6,*) ' Enter the coordinates of the field point'
      write (6,*) '                               99 to quit'
      write (6,*) ' ----------------------------------------'
      read  (5,*) x,y

      If(x.eq.99) Go to 99
      If(y.eq.99) Go to 99

      call lgf_2d_fs
     +              (Iopt
     +              ,x,y
     +              ,x0,y0
     +              ,G
     +              ,Gx,Gy
     +              )

      write (6,*) ' ---------------------------'
      write (6,*) ' Green function and gradient'
      write (6,*)
      write (6,100) G
      write (6,*)
      write (6,100) Gx,Gy
      write (6,*)
      write (6,*) ' ---------------------------'

c-------------------------------
c Test of integral identity
c by integrating: gradG.dot.normal_vector
c around a circle
c
c should be 0 or -1.0
c-------------------------------

      Else If(menu.eq.2) then

      write (6,*)
      write (6,*) " Enter: "
      write (6,*)
      write (6,*) "   the center of the test circle "
      write (6,*) "   the radius of the test circle"
      write (6,*) "   the number of integration points "
      write (6,*) " -----------------------------------"
      read (5,*) xcnt,ycnt,rad,mint

      mint1 = mint+1

      dth = pi2/(mint1-1.0D0)

      sm1 = 0.0D0

      Do i=1,mint

       th = (i-1.0D0)*dth
       cs = cos(th)
       sn = sin(th)

       x  = xcnt + rad*cs
       y  = ycnt + rad*sn

       call lgf_2d_fs
     +               (Iopt
     +               ,x,y
     +               ,x0,y0
     +               ,G
     +               ,Gx,Gy
     +               )

       vnx = cs    ! normal vector pointing outward from the circle
       vny = sn

       sm1 = sm1 + Gx*vnx + Gy*vny

      End Do

      sm1 = sm1 * dth * rad

      write (6,*)
      write (6,*) "Should be 0 or -1"
      write (6,*)
      write (6,100) sm1
      write (6,*)

c-----------
      End If
c-----------

      Go to 98

c-----
c Done
c-----

  99  Continue
 
 100  Format (3(2x,f15.10))
 110  Format (1x,i3,1x,f15.10,1x,f20.5)

      Stop
      End
