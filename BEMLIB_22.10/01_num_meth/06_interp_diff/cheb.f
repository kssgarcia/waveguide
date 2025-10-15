      program Chebyshev

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c        C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c----------------------------------------
c Applications of Chebyshev interpolation
c and approximation
c----------------------------------------

      Implicit Double Precision (a-h, o-z)

      Dimension x(1000),y(1000),c(0:200),d(0:200)
      Dimension Cheb(0:60,0:60)

c----------
c constants
c----------

      pi  = 3.14159265358D0
      pi2 = 2.0D0*pi

      Null = 0

c--------
c prepare
c--------

      open (7,file="PLOTDAT")  ! output file

 97   Continue

      write (6,*)
      write (6,*) " Menu of interpolated functions"
      write (6,*) 
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 for f(x) = x  "
      write (6,*) "       2 for f(x) = x**2 "
      write (6,*) "       3 for f(x) = exp(2x) "
      write (6,*) " ----------------------------"
      read  (5,*) Menu

      If(Menu.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Enter the polynomial order N "
      write (6,*) " -----------------------------"
      read  (5,*) N

      write (6,*)
      write (6,*) " Enter the end-points a, b"
      write (6,*) " -------------------------"
      read  (5,*) a,b

      Icomp = 2

c----
c prepare
c---

      N1 = N+1

c------------------------------
c Define data points
c
c abcsissas are zeros of T(N+1)
c------------------------------

      write (7,101) N1

      Do J=1,N1

       x(j) = Dcos((j-0.50D0)*pi/N1)

       If(Menu.eq.1) then
         y(j) = x(j)
       Else If(Menu.eq.2) then
         y(j) = x(j)**2
       Else If(Menu.eq.3) then
         y(j) = exp(2.0*x(j))
       End If

       write (7,101) j,x(j),y(j),y(j)

      End Do
   
c--------------------------
c Evaluate Cheb(xj)
c using a recursion relation
c--------------------------

      Do j=1,N1

        Cheb(0,j) = 1.0D0
        Cheb(1,j) = x(j)

        Do i = 2,N1
          Cheb(i,j) = 2.0D0*x(j)*Cheb(i-1,j)-Cheb(i-2,j)
        End Do

c       Do i = 0,N1
c         write (6,101) i,(Cheb(i,j),j=1,N1)
c       End Do

      End Do

c--------------------------
c Evaluate coefficients
c using Clenshaw
c--------------------------

      Do i=0,N
       c(i) = 0.0D0
       Do j=1,N1
        c(i) = c(i)+y(j)*Cheb(i,j)
       End Do
       c(i) = 2.0*c(i)/N1
      End Do

      c(0) = 0.5D0*c(0)

c----------
c Verify interpolation
c Clenshaw algorithm (6.4.15)
c----------

      write (6,*)
      write (6,*) " Specified and interpolated values"
      write (6,*)

      Do j=1,N1

        xx = x(j)
        d(N)   = c(N)
        d(N-1) = 2.0D0*xx*c(N)+c(N-1)
        Do i=N-2,0,-1
         d(i) = 2.0D0*xx*d(i+1)-d(i+2)+c(i)
        End Do
        value = d(0)-xx*d(1)
        write (6,101) j,xx,y(j),value
      End Do

c----------
c  Plotting
c----------

      write (6,*)
      write (6,*) " Please enter truncation level M "
      write (6,*)
      write (6,*) "         M=N gives interpolation "
      write (6,*) " --------------------------------"
      read  (5,*) M

      write (6,*)
      write (6,*) " How many plotting points ?"
      write (6,*) " --------------------------"
      read  (5,*) nprof

      dx = (b-a)/nprof
      nprof1 = nprof+1

      xx = a

      write (7,101) nprof1

c---
c Clenshaw
c---

      Do j=1,nprof1

        d(M)    = c(M)
        d(M-1) = 2.0*xx*c(M) + c(M-1)

        Do i=M-2,0,-1
         d(i) = 2.0D0*xx*d(i+1)-d(i+2)+c(i)
        End Do
        value = d(0)-xx*d(1)

        If(Menu.eq.1) then
         exact = xx
        Else If(Menu.eq.2) then
         exact = xx**2
        Else If(Menu.eq.3) then
         exact = exp(2.0*xx)
        End If

        write (6,101) j,xx,exact,value
        write (7,101) j,xx,exact,value

        xx = xx+dx

      End Do

      Go to 97  ! return to repeat

c-----
c Done
c-----

 99   Continue

      write (7,101) Null
      close (7)

 100  Format (30(1x,f10.5))
 101  Format (1x,I3,10(1x,f10.5))

      Stop
      End
