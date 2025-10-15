      program bernstein

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
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c-----------------------------------------
c Approximates a function with 
c Bernstein polynomials
c
c This program generates an array of values 
c of the approximated function
c for plotting purposes
c-----------------------------------------

      Implicit Double Precision (a-h, o-z)

      Dimension x(1000),y(1000)
      Dimension Bern(0:60,0:60)

      Parameter (k=200)

c----------
c constants
c----------

      pi  = 3.14159 265358
      pi2 = 2.0*pi

      null = 0

c----------
c launching
c----------

      write (6,*)
      write (6,*) "     MENU OF FUNCTIONS"
      write (6,*)
      write (6,*) " Enter 0 to quit"
      write (6,*) "       1 for f(x) = x  "
      write (6,*) "       2 for f(x) = x**2 "
      write (6,*) "       3 for f(x) = exp(2x) "
      write (6,*) " ----------------------------"

      read (5,*) Menu

      If(Menu.eq.0) Go to 99

      write (6,*) 
      write (6,*) " Data points will be spaced evenly"
      write (6,*) " over the interval (0, 1)"
      write (6,*) 
      write (6,*) " Enter number of intervals  M"
      write (6,*) " ----------------------------"
      read (5,*) M

      write (6,*) 
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) " 0 to quit "
      write (6,*) " 1 to use the recursive formula"
      write (6,*) "   of evaluating the Bernstein polynomials"
      write (6,*) " 2 to use the raw form"
      write (6,*) " -----------------------------------------"
      read (5,*) Icomp

      If(Icomp.eq.0) Go to 99

c--------
c prepare
c--------

      open (7,file="PLOTDAT")

      a = 0.0    ! left end
      b = 1.0    ! right end

      M1 = M+1
      h  = (b-a)/(M1-1.0)

c------------------------
c Compute the data points
c------------------------

      Do i=1,M1

       x(i) = a +(i-1.0)*h

       If(Menu.eq.1) then
         y(i) = x(i)
       Else If(Menu.eq.2) then
         y(i) = x(i)**2
       Else If(Menu.eq.3) then
         y(i) = exp(2.0*x(i))
       End If

c      write (7,101) i,x(i),y(i),y(i)

      End Do
   
c--------------------
c  Will plot k+1 points
c--------------------

      k1 = k+1
      write (7,101) k1

      t = a
      dt = (b-a)/k

      Do 10 j = 1,k1

c-----
c Computation of the Bernstein polynomials
c using the recursion formula
c given in equations (8.2.8)
c----

      If(Icomp.eq.1) then

       Do l=0,M
         Do i=0,M+1
          Bern(l,i) = 0.0D0
         End Do
       End Do

       Bern(0,1) = 1.0D0

       Do l=1,M
        Do i=1,M+1
          Bern(l,i) = (1.0-t)*Bern(l-1,i)+t*Bern(l-1,i-1)
        End Do
c      write (6,101) l,(Bern(l,kk),kk=1,M1)
       End Do

c-----
c Computation of the Bernstein polynomials
c with factorials
c not recommended
c-----

       Else If(Icomp.eq.2) then

        Do i=1,M1

         fcm = 1.0
         Do k1 = 1,M
           fcm = fcm*k1
         End Do

         ia = i-1
         fcia = 1.0
         Do k1 = 1,ia
           fcia = fcia*k1
         End Do

         mia = M+1-i
         fcmia = 1.0
         Do k1 = 1,mia
           fcmia = fcmia*k1
         End Do

         Bern(M,i) = fcM/(fcia*fcmia)*t**ia * (1.0-t)**mia

c        write (6,101) i,Bern(M,i)

        End Do

       End If

c--------------------------------------
c Generate the approximating polynomial
c--------------------------------------

       appr = 0.0

       Do i = 1,M1
        appr = appr + y(i)*Bern(M,i)
       End Do

       If(Menu.eq.1) then
         exact = t 
       Else If(Menu.eq.2) then
         exact = t**2
       Else If(Menu.eq.3) then
         exact = exp(2.0*t)
       End If

       write (6,101) j,t,exact,appr
       write (7,101) j,t,exact,appr

       t = t+dt

  10  Continue

c-----
c Done
c-----

   99 Continue

      write (7,101) null
      close  (7)

 100  Format (1x,f10.5,1x,f10.5,1x,f10.5)
 101  Format (1x,I3,10(1x,f10.5))

      Stop
      End
