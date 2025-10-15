      program bspl_un

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c              Oxford University Press
c------------------------------------------------

c--------------------------------
c Uniform B-spline approximation
c
c Will compute and plot uniform B-splines
c
c Section 8.4 (p.403)
c
c SYMBOLS:
c --------
c
c k: number of plotting points between knots
c--------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (k=32)

c---
c constants
c---

      null = 0

c---
c launching
c---

 98   Continue

      write (6,*)
      write (6,*) " Enter 0  to quit"
      write (6,*) "       1  to draw a linear B-spline"
      write (6,*) "       2  to draw a quadratic B-spline"
      write (6,*) "       3  to draw a cubic B-spline"
      write (6,*) " -------------------------------------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

c---
c prepare
c---

      open (9,file="PLOTDAT")
      k1 = k+1
      h  = 1.0/(k1-1.0)      ! step size


c-----------------------
      If(menu.eq.1) then   ! linear spline function
c-----------------------

        ktot = 2*k+1   ! number of printed points

        write (6,100) ktot
        write (9,100) ktot

        ic = 1
        x  = -1.0

        Do i=1,k
          y = 1.0+x
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          ic = ic+1
          x = x + h
        End Do

        Do i=1,k1
          y = 1-x
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          x = x + h
        End Do

c----------------------------
      Else If(menu.eq.2) then      !   equations (8.4.8)
c----------------------------

        ktot = 3*k+1    ! number of printed points

        write (6,100) ktot
        write (9,100) ktot

        ic = 1
        x  = -1.0

        Do i=1,k
          y = 0.5*(1.0+x)**2
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          ic = ic+1
          x  = x + h
        End Do

        Do i=1,k
          y = 0.5*(1.0+2.0*x-2.0*x**2)
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          ic = ic+1
          x  = x + h
        End Do

        Do i=1,k1
          y = 0.5*(2.0-x)**2
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          ic = ic+1
          x  = x + h
        End Do

c----------------------------
      Else If(menu.eq.3) then     !  equations (8.4.14)
c----------------------------

        ktot = 4*k+1     ! number of printed points

        write (6,100) ktot
        write (9,100) ktot

        ic = 1
        x  = -2.0

        Do i=1,k
          y = (x+2.0)**3/6.0
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          ic = ic+1
          x  = x + h
        End Do

        Do i=1,k
          y = (4.0-6.0*x**2-3.0*x**3)/6.0
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          ic = ic+1
          x  = x + h
        End Do

        Do i=1,k
          y = (4.0-6.0*x**2+3.0*x**3)/6.0
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          ic = ic+1
          x  = x + h
        End Do

        Do i=1,k1
          y = (2.0-x)**3/6.0
          write (6,100) ic,x,y
          write (9,100) ic,x,y
          ic = ic+1
          x  = x + h
        End Do

c-----------
      End If
c-----------

c---
c repeat
c---

      Go to 98

c-----
c Done
c-----

  99  Continue

      write (9,100) null
      close (9)

 100  Format (1x,i3,1x,f15.8,1x,f15.8)

      Stop
      End
