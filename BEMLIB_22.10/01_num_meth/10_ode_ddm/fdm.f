      program fdm

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c   This program accompanies the book:
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c------------------------------------------------

c-------------------------------------
c Finite-difference solution of an ODE
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(200),b(200),c(200),s(200),x(200)

c---
c constants
c---

      pi = 3.14159265358
      pih = 0.50*pi

      null = 0
      none = 1

      zero =  0.0
      one  =  1.0
      onem = -1.0

c---
c prepare
c---
      open (1,file="PLOTDAT")

c---
c preferences
c---

  98  Continue

      write (6,*)
      write (6,*) " Please enter epsilon"
      write (6,*) "           99 to quit"
      write (6,*) " --------------------"
      read  (5,*) eps

      If(eps.eq.99) Go to 99

      write (6,*)
      write (6,*) " Please enter the number of intervals N"
      write (6,*) "                              0 to quit"
      write (6,*) " --------------------------------------"
      read  (5,*) n

      If(n.eq.0) go to 99

c---
c preparations
c---

      n1 = n+1
      h  = pih/(n1-1.0)
      hs = h**2

c---
c tridiagonal matrix (10.3.9)
c---

      Do i=1,n
       b(i) = 1.0
       c(i) = 1.0
       ti   = (i-1.0)*h
       a(i) = -2.0+hs*(4.0+eps*cos(ti))
       s(i) = 0.
      End Do

      b(1) = 2.0

c---
c right-hand side, eq (10.3.10)
c---

      s(1) = -4.0*h
      s(n) = 1.0

c---
c solve the system
c---

      call thomas (n,a,b,c,s,x)

c---
c  print
c---

      write (1,*) n1

      If(eps.eq.0) then
         write (6,*)
         write (6,*) " t, x, exact, error"
         write (6,*)
      Else
         write (6,*)
         write (6,*) " t, x"
         write (6,*)
      End If


      Do i=1,n

       ti = (i-1.0D0)*h

       If(eps.eq.0) then
         ti2 = 2.0*ti
         exact = cos(ti2)-sin(ti2)   ! exact solution for eps = 0
         error = x(i)-exact
         write (6,100) i,ti,x(i),exact,error
         write (1,100) i,ti,x(i),exact,error
       Else 
         write (6,100) i,ti,x(i)
         write (1,100) i,ti,x(i)
       End If

      End Do

      If(eps.eq.0) then
        error = 0.0
        write (6,100) n1,pih,onem,onem,error
        write (1,100) n1,pih,onem,onem,error
      Else
        write (6,100) n1,pih,onem
        write (1,100) n1,pih,onem
      End If

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
