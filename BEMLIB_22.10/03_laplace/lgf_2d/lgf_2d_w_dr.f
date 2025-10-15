      program lgf_2d_w_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------
c Driver for the Green's function of Laplace's
c equation in a domain bounded by a plane wall
c located at y = wall
c-------------------------------------------
 
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

      Iopt = 2   ! will need G and grad(G)

c-------
c launch
c-------

      write (6,*)
      write (6,*) '         MENU'
      write (6,*)
      write (6,*) ' Enter 0 to quit'
      write (6,*) '       1 for one evaluation'
      write (6,*) '       2 to test integral identities'
      write (6,*) '------------------------------------'
      read  (5,*) menu

      If(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) ' The wall is located at y = wall'
      write (6,*)
      write (6,*) ' Please enter wall'
      write (6,*) '                  99 to quit'
      write (6,*) ' ---------------------------'
      read  (5,*) wall

      If(wall.eq.99) Go to 99

      write (6,*) ' Enter 0 to quit'
      write (6,*) '       1 for the Green function'
      write (6,*) '       2 for the Neumann function'
      write (6,*) '---------------------------------'
      read  (5,*) Ign

  98  Continue

      write (6,*)
      write (6,*) ' Enter the x, y coordinates of one singularity'
      write (6,*) '           99 for either to quit'
      write (6,*) ' -------------------------------'
      read  (5,*) x0,y0

      If(x0.eq.99) Go to 99
      If(y0.eq.99) Go to 99

c---------------
c One evaluation
c---------------

      IF(menu.eq.1) then

      write (6,*)
      write (6,*) ' Enter the (x,y) of the velocity point'
      write (6,*) '                            99 to quit'
      write (6,*) ' -------------------------------------'
      read  (5,*) x,y

      If(x.eq.99) Go to 99
      If(y.eq.99) Go to 99

      call lgf_2d_w 
     +              (Iopt
     +              ,Ign
     +              ,x,y
     +              ,x0,y0
     +              ,wall
     +              ,G
     +              ,Gx,Gy
     +              )

      write (6,*) ' ----------------------------'
      write (6,*) ' Green function and gradient:'
      write (6,100) G
      write (6,100) Gx,Gy
      write (6,*) ' ---------------------------'

c----------------------------
c Test of integral properties
c
c sm1 should be 0 or -1
c----------------------------

      Else If(menu.eq.2) then

      write (6,*)
      write (6,*) " Please enter:"
      write (6,*)
      write (6,*) "  the center of the test circle "
      write (6,*) "  the radius of  the test circle"
      write (6,*) "  the number of integration points "
      write (6,*) " ----------------------------------"
      read  (5,*) xcnt,ycnt,rad,mint

      dth = pi2/mint

      sm1 = 0.0D0

      Do  i=1,mint

        th = (i-1.0D0)*dth
        cs = cos(th)
        sn = sin(th)

        x  = xcnt + rad*cs
        y  = ycnt + rad*sn

       call lgf_2d_w 
     +              (Iopt
     +              ,Ign
     +              ,x,y
     +              ,x0,y0
     +              ,wall
     +              ,G
     +              ,Gx,Gy
     +              )

        vnx = cs    ! normal vector pointing outward
        vny = sn

        sm1 = sm1 + Gx*vnx + Gy*vny

      End Do

      cf  = dth*rad
      sm1 = sm1 * cf

      write (6,*)
      write (6,*) " Should be 0 or -1"
      write (6,100)
      write (6,100) sm1
      write (6,100)

c--------
      End If
c--------

      Go to 98

c-----
c Done
c-----

  99  Continue
 
 100  Format (3(2x,f15.10))
 110  Format (1x,i3,1x,f15.10,1x,f20.5)

      Stop
      End
