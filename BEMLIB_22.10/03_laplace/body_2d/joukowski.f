      program joukowski

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------
c Generate Joukowski airfoils
c
c according to equation (7.10.2) of Pozrikidis (1997, p.372)
c
c lamda will be set to 1.0
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c constants
c---

      pi = 3.14159265358

      null= 0

c---
c parameters
c---

      write (6,*) 
      write (6,*) "  MENU"
      write (6,*) 
      write (6,*) "  Enter "
      write (6,*) 
      write (6,*) "  0 to quit "
      write (6,*) "  1 for a Joukowski airfoil"
      write (6,*) "  ---------------------------------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

      open (2,file="joukowski.xy")

c----------------
      If(Menu.eq.1) then
c----------------

  98  Continue

      write (6,*)
      write (6,*) "  Please enter the number of points"
      write (6,*) "  0 to quit"
      write (6,*) "  ---------"
      read  (5,*) n

      If(n.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter x and y coordinates of circle center"
      write (6,*) " 99 to quit"
      write (6,*) " ----------"
      read (5,*) xc,yc

      If(xc.eq.99.or.yc.eq.99) Go to 99

c---
c compute the airfoil
c---

      n1     = n+1
      xstart = -pi
      xend   =  pi
      step   = (xend-xstart)/(n1-1.0)
      Rad2   = (xc+1.0)**2+yc**2
      Rad    = sqrt(Rad2)

      th = xstart

      write (2,900) n1,xc,yc
      write (6,900) n1,xc,yc

      Do i=1,n1
        x = xc+Rad*cos(th)
        y = yc+Rad*sin(th)
        r2 = x**2+y**2
        xx = x + x/r2
        yy = y - y/r2
        write (2,900) i,xx,yy
        write (6,900) i,xx,yy
        th = th +step
      End Do

      Go to 98

c------------
      End If
c------------

  99  Continue

      write (2,900) null
      close (2)

c-----
c Done
c-----

 900  Format (1x,i3,4(1x,f10.5))

      Stop
      End
