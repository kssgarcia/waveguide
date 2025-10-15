      program thomas_dr

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
c           C. Pozrikidis
c
c ``Numerical Computation in Science and Engineering''
c
c        Oxford University Press, 1998
c------------------------------------------------

c------------------------------
c Driver for Thomas' algorithm
c for tridiagonal systems
c
c Algorithm (3.4.1)
c-------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(1026),b(1026),c(1026),s(1026),x(1026)

c----------
c constants
c----------

      pi = 3.14159 265358 D0

c------------
c preferences
c------------

      write (6,*)
      write (6,*) " Choose problem to be solved"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " 1 to read data from file thomas.dat"
      write (6,*) "                       or thomas_pr.dat"
      write (6,*) " 2 to solve the ode described in text   "
      write (6,*) " ---------------------------------- "
      read  (5,*) menu

      If(menu.eq.0) Go to 99

 98   Continue

c---
      If(menu.eq.1) then
c---

       write (6,*) " Enter:"
       write (6,*)
       write (6,*) " 0 to quit"
       write (6,*) " 1 for the standard algorithm "
       write (6,*) " 2 for the compact algorithm "
       write (6,*) " 3 for periodic boundary conditions"
       write (6,*) " ---------------------------------- "
       read  (5,*) method

c---
      Else
c---

      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 to quit"
      write (6,*) " 1 for the standard algorithm "
      write (6,*) " 2 for the compact algorithm "
      write (6,*) " ---------------------------- "
      read  (5,*) method

c---
      End If
c---

      If(method.eq.0) Go to 99

c-----------------------
      If(menu.eq.1) then
c-----------------------

      If(method.eq.1.or.method.eq.2) then

        open (8,file="thomas.dat")
        read (8,*) n
        read (8,*) (a(i),i=1,n)
        read (8,*) (b(i),i=1,n-1)
        read (8,*) (c(i),i=2,n)
        read (8,*)
        read (8,*) (s(i),i=1,n)
        close(8)

      Else If(method.eq.3) then          ! periodic BC

        open (8,file="thomas_pr.dat")
        read (8,*) n
        read (8,*) (a(i),i=1,n)
        read (8,*) (b(i),i=1,n)
        read (8,*) (c(i),i=1,n)
        read (8,*)
        read (8,*) (s(i),i=1,n)
        close(8)

      End If

c----
c call thomas
c----

      If(method.eq.1) then
          call thomas (N,a,b,c,s,x)

      Else If(method.eq.2) then
          call thomas_c  (N,a,b,c,s)

      Else If(method.eq.3) then
          call thomas_pr (N,a,b,c,s,x)

      End If

c----
c Display
c----

      write (6,*)
      write (6,*) "Solution vector:"
      write (6,*)

      Do i=1,N
         If(menu.eq.1) write (6,100) i,x(i)
         If(menu.eq.2) write (6,100) i,s(i)
      End Do

c----------------------------
      Else If(menu.eq.2) then
c----------------------------

      write (6,*) 
      write (6,*) " Enter the length of the solution domain"
      write (6,*) "       in multiples of pi"
      write (6,*) 
      write (6,*) "       0 to quit"
      write (6,*) " ---------------"
      read (5,*) bb

      If(bb.lt.0.0000000001) Go to 99

      write (6,*)
      write (6,*) " Enter number of intervals: N"
      write (6,*) 
      write (6,*) "       0 to quit"
      write (6,*) " ---------------"
      read  (5,*) N

      If(N.eq.0) Go to 99

      Na = N-1
      N1 = N+1

      aa = 0.0
      bb = bb*pi
      dx = (bb-aa)/N
      cn = 2.0-4.0*dx**2

      Do i=1,Na
        a(i) =  cn
        b(i) = -1.0D0
        c(i) = -1.0D0
        s(i) =  0.0D0
      End Do

      s(1)  = 0
      s(Na) = 5.0D0

c----
c call thomas
c---

      If(method.eq.1) call thomas (Na,a,b,c,s,x)
      If(method.eq.2) call thomas_c (Na,a,b,c,s)

c----
c Display
c---

      x(N1) = 5.0

      Do i=Na,1,-1
        If(method.eq.1) x(i+1) = x(i)
        If(method.eq.2) x(i+1) = s(i)
      End Do

      x(1) = 0.0

      write (6,*)
      write (6,*) "Solution vector:"
      write (6,*)

      Do i=1,N1
        write (6,100) i,x(i)
      End Do

c-----------
      End If   ! loop over menu
c-----------

      Go to 98   ! repeat

 99   Continue

      write (6,*)
      write (6,*) "Thank you for your time"

c-----
c Done
c-----

 100  Format (1x,i3,1x,f15.8)

      Stop
      End
