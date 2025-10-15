      program lagrange_int

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
c            C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press, 1998
c------------------------------------------------

c---------------------------------------
c Applications of Lagrange interpolation
c---------------------------------------

      Implicit Double Precision (a-h, o-z)

      Dimension x(-1000:1000),f(-1000:1000)

      Parameter (nprof=100)    ! number of plotting points

      Dimension Z(20)

c----------
c constants
c----------

      pi  = 3.141592 65358 D0
      pi2 = 2.0D0*pi

      null = 0
      zero = 0.0

c--------
c prepare
c--------

      open (7,file="PLOTDAT")

      nprof1 = nprof+1

c------------
c preferences
c------------

 98   Continue

      write (6,*) 
      write (6,*) "            MENU"
      write (6,*) 
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 to generate the function 2**n Phi"
      write (6,*) "   with evenly spaced points in [-1, 1]"
      write (6,*) " 2 for the Lagrange polynomials"
      write (6,*) " 3 for Lagrange interpolation "
      write (6,*) "------------------------------"
      read  (5,*) menu

      If(menu.eq.0) Go to 99

      write (6,*)
      write (6,*) " Enter:"
      write (6,*) 
      write (6,*) " 0 to quit"
      write (6,*) " 1 to generate data at evenly-spaced points" 
      write (6,*) " 2 to generate data at Lobatto base-points" 
      write (6,*) " 3 to generate data at Cheb-sk base-points" 
      write (6,*) " 9 to read from file: int_1d.dat" 
      write (6,*) "-------------------------------------"
      read  (5,*) Idata

c------------------------
      If(Idata.eq.9) then       ! read the data
c------------------------

       open(3,file="int_1d.dat")

         read (3,*) N1
         Do i=1,N1
           read (3,*) idum,x(i),f(i)
         End Do

         N = N1-1

       close (3)

       Go to 93

c-----------
      End If
c-----------

      write (6,*)
      write (6,*) " Enter first and last points: a, b"
      write (6,*) " ---------------------------------"
      read  (5,*) a,b

      write (6,*)
      write (6,*) " Enter the polynomial order, N"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read (5,*) N 

      If(N.eq.0) Go to 99

      N1 = N+1

c--------------------
       If(Idata.eq.1) then   ! evenly-spaced points
c---------------------

          h = (b-a)/(N1-1.0D0)
          Do i=1,N1
            x(i) = a +(i-1.0D0)*h
          End Do

c------------------
        Else If(Idata.eq.2) then   ! Lobatto points
c-------------------

          call lobatto (N,Z)

          Do i=1,N1
            x(i) = (a+b)/2.0D0 + Z(i)*(b-a)/2.0D0
          End Do

c--------------------
       Else If(Idata.eq.3) then   ! second-kind Cheb second kind
c---------------------

          h = pi/(N1-1.0D0)
          Do i=1,N1
            Z(i) = Dcos((i-1.0D0)*h)
            x(i) = (a+b)/2.0D0 + Z(i)*(b-a)/2.0D0
          End Do

c-------------
        End If
c-------------

        If(menu.eq.3) then

        write (6,*)
        write (6,*) " Choose the function"
        write (6,*)
        write (6,*) " Enter 0 to quit"
        write (6,*) "       1 for f(x) = 1/(1+25*x**2) "
        write (6,*) "       2 for f(x) = exp(2x) "
        write (6,*) "       3 for f(x) = exp(3x)/(1+25*x**2)"
        write (6,*) "       4 for f(x) = sqrt(x)"
        write (6,*) "       5 for f(x) = 1/(1+4*x**2) "
        write (6,*) "  -------------------------------------"
        read (5,*) Ifunc

        If(Ifunc.eq.0) Go to 99

        Do i=1,n1
          If(ifunc.eq.1) then
             f(i) = 1.0D0/(1.0D0+25.0D0*x(i)**2)
          Else If (Ifunc.eq.2) then
             f(i) = exp(2.0D0*x(i))
          Else If (Ifunc.eq.3) then
             f(i) = exp(3.0D0*x(i))/(1.0+25.0D0*x(i)**2)
          Else If (Ifunc.eq.4) then
             f(i) = Dsqrt(x(i))
          Else If (Ifunc.eq.5) then
             f(i) = 1.0D0/(1.0D0+4.0D0*x(i)**2)
          End If
        End Do

        End If

c------------
c END OF DATA
c------------

  93  Continue

c-----------------
c Display the data
c-----------------

      write (6,*)
      write (6,*) " Data points:"
      write (6,*) " -----------"

      Do i=1,n1
       write (6,101) i,x(i),f(i)
      End Do

c---------------
c Prepare to run
c---------------

      If(menu.eq.2) then    ! will generate the Lagrange polynomials

        write (6,*)
        write (6,*) " Please enter the number of the polynomial"
        write (6,*) " -----------------------------------------"
        read (5,*) l1

      End If

c----------
c  Plotting
c----------

      xx = a
      dx = (b-a)/nprof

      write (7,101) nprof1

      Do j=1,nprof1
       
c-----------------------
      If(menu.eq.1) then    ! will produce the generating function
c-----------------------

        prod = 1.0
        Do l2=1,n1
         prod = prod*(xx-x(l2))
        End Do

        sum = 2.0**n * prod

        write (7,101) j,xx,sum

c-----------------------------
      Else If (menu.eq.2) then  ! will produce the Lagrange polynomials
c-----------------------------

        prod = 1.0D0

        Do l2=1,n1
          If(l2.ne.l1) prod = prod*(xx-x(l2))/(x(l1)-x(l2))
        End Do

        sum = prod

        arg = n*pi*(xx-x(l1))/(b-a)
        sonc = sin(arg)/arg

        write (7,101) j,xx,sum,sonc

c-----------------------------
      Else if (menu.eq.3) then  ! Lagrange interpolation
c-----------------------------

        sum = 0.0

        Do l1=1,n1
          prod = 1.0
          Do l2=1,n1
            If(l2.ne.l1) prod = prod*(xx-x(l2))/(x(l1)-x(l2))
          End Do
         sum = sum + prod*f(l1)
        End Do

c---
c exact value
c---

        If(Ifunc.eq.1) then
         exact = 1.0D0/(1.0D0+25.0D0*xx**2)
        Else If(Ifunc.eq.2) then
         exact = exp(2.0*xx)
        Else If(Ifunc.eq.3) then
         exact = exp(3.0*xx)/(1.0+25.0*xx**2)
        Else If(Ifunc.eq.4) then
         exact = sqrt(xx)
        Else If(Ifunc.eq.5) then
         exact = 1.0D0/(1.0D0+4.0D0*xx**2)
        End If

        error = sum-exact

        write (7,101) j,xx,exact,sum,error

c-----------
      End If
c-----------

      xx = xx+dx

      End Do

      Go to 98       ! repeat

c-----
c Done
c-----

 99   Continue

      write (7,101) null
      Close (7)

 100  Format (1x,f10.5,1x,f10.7,1x,f10.7)
 101  Format (1x,I3,10(1x,f10.7))

      Stop
      End
