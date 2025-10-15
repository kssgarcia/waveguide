      program fft_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c              C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press
c------------------------------------------------

c-------------------------------------------
c Driver for the Fast Fourier Transform: fft
c-------------------------------------------

      Implicit Double Precision (a-h, o-z)

      Dimension x(4096),f(4096,2),z(4096)
      Dimension g(0:4096,2)
      Dimension af(0:4096),bf(0:4096),cfm(0:4096)

c----------
c constants
c----------

      pi = 3.14159 265358 D0
      pi2 = 2.0D0*pi

      null = 0

c------------
c preferences
c------------

  98  Continue

      write (6,*)
      write (6,*) "    MENU OF FUNCTIONS"
      write (6,*)
      write (6,*) " Enter 1 for f(x) = 1.0  "
      write (6,*) "       2 for f(x) = x**2 "
      write (6,*) "       3 for f(x) = exp(x) "
      write (6,*) "       4 for f(x) = exp(cos(2*pi*x) "
      write (6,*) "       5 for f(x) = exp(cos(  pi*x) "
      write (6,*) "       6 for f(x) = cos(2*pi*x) "
      write (6,*) "       7 for f(x) = sin(2*pi*x) "
      write (6,*) "     100 to read the data from file: fft.inp"
      write (6,*) "       0 to quit  "
      write (6,*) " -----------------"
      read (5,*) menu

      if(menu.eq.0) Go to 99

c-------------------------
c generate the data points
c-------------------------

      if(menu.lt.100) then

       write (6,*) 
       write (6,*) " Enter limits of interpolation a, b"
       write (6,*) " ----------------------------------"
       read  (5,*) a,b
c      a=0.0D0
c      b=1.0D0

       write (6,*) 
       write (6,*) " Enter exponent m"
       write (6,*) 
       write (6,*) " The number of data points is: n = 2**m"
       write (6,*) " --------------------------------------"
       read  (5,*) nexp
c      nexp=3

       n = 2**nexp

       write (6,*)
       write (6,*) " Number of data points : ",n

       n1 = n+1
       h = (b-a)/(n1-1.0D0)
 
       Do i=1,n
        x(i) = a +(i-1.0D0)*h
        if(menu.eq.1) f(i,1) = 1.0
        if(menu.eq.2) f(i,1) = x(i)**2
        if(menu.eq.3) f(i,1) = exp(x(i))
        if(menu.eq.4) f(i,1) = exp(cos(pi2*x(i)))
        if(menu.eq.5) f(i,1) = exp(cos(pi *x(i)))
        if(menu.eq.6) f(i,1) = cos(pi2*x(i))
        if(menu.eq.7) f(i,1) = sin(pi2*x(i))
        f(i,2) = 0.0D0
       End Do

      End If

c-----------------
c read data points
c-----------------

      If(menu.eq.100) then

c      open (2,file="fft.dat")
       open (2,file="keller.dat")

       read (2,*) nexp,h

       n = 2**nexp

       write (6,*)
       write (6,*) " fft_dr: number of data points: ",n

       Do i=1,n
         x(i) = (i-1.0D0)*h
         read (2,*) idle,f(i,1)
         f(i,2) = 0.0D0          ! complex part set to zero
       End Do

       close (2)

      End If

c-------------------------------------
c reduce according to equation (8.7.1)
c-------------------------------------

      n1 = n+1
      fc = 1/(n1-1.0D0)

      Do i=1,n
       g(i-1,1) = fc*f(i,1)
       g(i-1,2) = fc*f(i,2)
      End Do

c-----------
c Do the fft
c-----------

      call fft 
     +
     +  (nexp
     +  ,g
     +  ,af,bf
     +  )

c---------
c printing
c---------

      open (1,file="fft.out")

      write (6,*)
      write (6,*) " Fourier Coefficients"
      write (6,*)
      write (6,*) "  i, a(i), b(i), |c(i)|"
      write (6,*) " ----------------------"
      write (6,*)

      write (1,101) n
      write (6,101) n

      Do i=0,n-1

        cfm(i) = sqrt(af(i)**2+bf(i)**2)   ! magnitude of power spectrum

        write (6,101) i,af(i),bf(i),cfm(i)
        write (1,101) i,af(i),bf(i),cfm(i)
      End Do

      write (1,101) null

      close (1)

c----------------------
c Confirm interpolation
c----------------------

      write (6,*)
      write (6,*) " Verified values"
      write (6,*)
      write (6,*) " i, x, f, F"
      write (6,*) " ----------"

      nh = n/2
      wn = pi2/(n*h)   ! fundamental wave number

      Do i=1,n

        xh = x(i)-x(1)     ! hat{x}

        z(i) = 0.5D0*af(0)

        Do j=1,n/2
         arg  = j*wn*xh
         cs   = Dcos(arg)
         sn   = Dsin(arg)
         z(i) = z(i)+af(j)*cs+bf(j)*sn
        End Do

        write (6,101) i,x(i),f(i,1),z(i)

      End Do

c-------
c repeat
c-------

      Go to 98

c-----
c done
c-----

  99  Continue

 100  Format (1x,f10.5,1x,f10.5,1x,f10.5)
 101  Format (1x,I8,10(1x,f10.5))

      Stop
      End
