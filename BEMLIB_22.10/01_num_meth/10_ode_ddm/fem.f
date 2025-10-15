      program fem 

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c----------------------------------------
c Finite-element solution of equation (10.4.1)
c with boundary conditions (10.4.2)
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(200),b(200),c(200),s(200),x(200)

c---
c constants
c---

      null = 0
      none = 1

      zero = 0.0
      one  = 1.0

c---
c prepare
c---

      open (1,file="fem.out")

c---
c input
c---

  98  Continue

      write (6,*)
      write (6,*) " Please enter the number of intervals N"
      write (6,*) "              0 to quit"
      write (6,*) " ----------------------"
      read  (5,*) n

      If(n.eq.0) Go to 99

c---
c prepare
c---

      n1 = n+1
      h  = 1.0/(n1-1.0)

c---
c linear system (10.4.10)
c---

      aa = 2.0/3.0 * (3.0/h-h)
      bb = -1/h-h/6.0

      Do i=1,n
       a(i) = aa
       b(i) = bb
       c(i) = bb
       ti   = i*h
       s(i) = h*ti
      End Do
      a(n) = 0.5 * a(n)
      s(n) = 0.5 * h*(1-h/3.0)

      call thomas (n,a,b,c,s,x)

c---
c print numerical and exact solution
c---

      write (6,*)
      write (6,*) " t, x, exact, error"
      write (6,*)
      write (1,*) n1

      write (6,100) none,zero,zero,zero
      write (1,100) none,zero,zero,zero

      Do i=1,n
       ti    = i*h
       exact = sin(ti)/cos(one)-ti
       error = x(i)-exact
       write (6,100) i+1,ti,x(i),exact,error
       write (1,100) i+1,ti,x(i),exact,error
      End Do

      Go to 98

c-----
c Done
c-----

 99   Continue

      write (1,*) null
      close (1)

 100  Format (1x,I3,10(1x,f12.7))
 101  Format (30(1x,f10.7))

      Stop
      End
