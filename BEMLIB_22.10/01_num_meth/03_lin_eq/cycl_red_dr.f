      program cycl_red_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press
c
c 1998
c------------------------------------------------

c---------------------------
c Driver for algorithm 3.4.2
c Cyclic reduction
c---------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension s(200),x(200)
      Integer p

c---
c constants
c---

      pi = 3.14159 265358

c---
c launch
c---

 98   Continue

      write (6,*)
      write (6,*) " Choose problem to be solved"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 to read data from a file"
      write (6,*) "       2 to solve ode described in text "
      write (6,*) "       -------------------------------- "
      read  (5,*) menu

      If(menu.eq.0) Go to 99

c-----------------------
      If(menu.eq.1) then
c-----------------------

        open (8,file="cycl_red.dat")

        read (8,*) p
        n = 2**p
        read (8,*) a,b,c
        read (8,*)
        read (8,*) (s(i),i=1,n)

        close(8)

c-----------------------
      Else If(menu.eq.2) then
c-----------------------

        write (6,*) 
        write (6,*) " Specify the length of the interval"
        write (6,*) " in multiples of pi"
        write (6,*) "                   0 to quit"
        write (6,*) " -----------------------------------"
        read  (5,*) bb

        If(bb.eq.0) Go to 99

        bb  = bb*pi
        aa  = 0
        b   = -1.0
        c   = -1.0

        write (6,*) 
        write (6,*) " Please enter the exponent p"
        write (6,*) "                   0 to quit"
        write (6,*) " ---------------------------"
        read  (5,*) p

        If(p.eq.0) Go to 99

        m   = 2**p+1
        dx  = (bb-aa)/m
        dx2 = dx**2
        n   = m-1
        a   = 2.0-4.0*dx2

        s(1) = 0
        Do i=2,n-1
         s(i) = 0
        End Do
        s(n) = 5.0

c-----------
      End If
c-----------

c---------------------------
c proceed with the reduction
c---------------------------

      call cycl_red (n,p,a,b,c,s,x)

c---------
c printing
c---------

      Do i=1,n
        write (6,100) i,x(i)
      End do
 
c----------
c verifying
c----------

      write (6,*)
      write (6,*) " Residuals"
      write (6,*)

      i = 1
      res = s(1)-a*x(1) -b*x(2)
      write (6,100) i,res

      Do i=2,n-1
        res = s(i)-c*x(i-1)-a*x(i)-b*x(i+1)
        write (6,100) i,res
      End Do

      res = s(n) -a*x(n)-c*x(n-1)
      write (6,100) i,res

      Go to 98   ! return to repeat

c-----
c Done
c-----

 99   Continue

      write (6,*) "Thank you for running me"

 100  Format (1x,i3,1x,f10.6)

      Stop
      End
