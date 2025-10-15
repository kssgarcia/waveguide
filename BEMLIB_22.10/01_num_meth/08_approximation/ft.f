      program ft

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c------------------------------------
c Complex, cosine, and sine transform
c of a function of one variable f(x)
c computed by the conventional method
c
c Points are spaced evenly over [a, b]
c
c SYMBOLS:
c -------
c
c a:  left end of interpolation domain
c b:  right end of interpolation domain
c L:  length of the domain of interpolation, L = b-a
c h:  point spacing
c
c af, bf:	real and imaginary parts of complex coefficients
c
c ac:    coefficinets of cosine Fourier transform
c as:	 coefficients of sine Fourier transform
c------------------------------------

      Implicit Double Precision (a-h, o-z)

      Double Precision L,k,kh

      Integer p

      Dimension x(32000),f(32000),z(32000),w(32000),v(32000)
      Dimension af(0:32000),bf(0:32000)
      Dimension ac(0:32000),as(0:32000)

c----------
c constants
c----------

      pi  = 3.14159 265358
      pi2 = 2.0*pi

      null = 0

c---
c parameters
c---

  98  Continue

      write (6,*)
      write (6,*) " MENU OF FUNCTIONS"
      write (6,*)
      write (6,*) " Enter:"
      write (6,*)
      write (6,*) "  1 for f(x) = 1.0  "
      write (6,*) "  2 for f(x) = x**2 "
      write (6,*) "  3 for f(x) = exp(x) "
      write (6,*) "  4 for f(x) = exp(cos(2*pi*x) "
      write (6,*) "  5 for f(x) = exp(sin(  pi*x) "
      write (6,*) " 10 to read the data from file: ft.dat"
      write (6,*) " 20 for some data"
      write (6,*) "  0 to quit  "
      write (6,*) "-------------"
      read (5,*) menu

      If(menu.eq.0) Go to 99

c--------------------
c generate data points
c--------------------

      if(menu.lt.10) then

       write (6,*)
       write (6,*) " Enter the end-points a, b"
       write (6,*) " -------------------------"
       read  (5,*) a,b

       write (6,*)
       write (6,*) " Enter the number of intervals N"
       write (6,*) " -------------------------------"
       read  (5,*) N

       N1 = N+1
       L  = b-a
       h  = L/(N1-1.0D0)

       Do i=1,N1
        x(i) = a+(i-1.0)*h
        if(menu.eq.1) f(i) = 1.0D0
        if(menu.eq.2) f(i) = x(i)**2
        if(menu.eq.3) f(i) = exp(x(i))
        if(menu.eq.4) f(i) = exp(cos(pi2*x(i)))
        if(menu.eq.5) f(i) = exp(sin( pi*x(i)))
       End Do

      End If

c----------
c read data
c-----------

      if(menu.eq.10) then

       open (2,file="ft.dat")

       read (2,*) N,h

       write (6,*)
       write (6,*) " Number of intervals: ",N

       N1 = N+1

       fmean = 0.0D0

       Do i=1,N1
         x(i) = (i-1.0D0)*h
         read (2,*) idle,a,a,f(i)
         fmean = fmean+f(i)
c        write (6,*) idle,a,a,f(i)
       End Do

       fmean = fmean/N1

       Do i=1,N1
        f(i) = f(i)-fmean
       End Do

       close (2)

       L = N*h
       write (6,*) " Period: ",L

      end if

c----------
c some data
c-----------

      if(menu.eq.20) then
          N = 8
          f(1) = -0.12096097061131
          f(2) =  0.03073269850607
          f(3) =  0.16442356964692
          f(4) =  0.20179734366240
          f(5) =  0.12096097061131
          f(6) = -0.03073269850607
          f(7) = -0.16442356964692
          f(8) =  -0.20179734366240
          f(9) =  -0.12096097061131

          f(1) =  -0.07548045967823
          f(2) =   0.08073424015518
          f(3) =   0.18965591705357
          f(4) =   0.18747972992629
          f(5) =   0.07548045967823
          f(6) =  -0.08073424015518
          f(7) =  -0.18965591705357
          f(8) =  -0.18747972992629
          f(9) =  -0.07548045967823

          N1 = N+1
          h = 1.0;
          L = N*h
      end if

c--------
c prepare
c--------

      open (7,file="ft.out")

      k  = pi2/L
      kh = k*h
      ff = 2.0D0/N
 
c-------------------------------
c Even/Odd number of intervals ?
c-------------------------------

      If(mod(N,2).eq.0) then
        write (6,*) N," is even"
        M  = N/2
        fc = 0.50D0
      Else
        write (6,*) N," is odd"
        M = (N-1)/2
        fc = 1.0D0
      End If

c--------
c prepare
c--------

      M2 = 2*M
      M4 = 4*M

      NFC = M4   ! number of fourier coefficients

      NFC = 20

c-----------------------------
c Compute the Fourier coefficients
c-----------------------------

      Do p=0,NFC         ! loop over coefficients

       af(p) = 0.5D0*(f(1)+f(N1))
       bf(p) = 0.0D0

       ac(p) = 0.5D0*(f(1)+f(N1)*cos(p*pi))
       as(p) = 0.0D0

       Do i=2,N

c---
c complex transform
c---

        arg   = p*(i-1.0D0)*kh
        cs    = Dcos(arg)
        sn    = Dsin(arg)
        af(p) = af(p) + cs*f(i)
        bf(p) = bf(p) + sn*f(i)

c---
c cosine and sine transforms
c---
        arg = 0.5D0*arg
        cs  = Dcos(arg)
        sn  = Dsin(arg)

        ac(p) = ac(p) + cs*f(i)
        as(p) = as(p) + sn*f(i)

       End Do

c---
c normalize
c---

       af(p) = ff*af(p)
       bf(p) = ff*bf(p)

       ac(p) = ff*ac(p)
       as(p) = ff*as(p)

      End Do

      af(M) =   fc*af(M)
      ac(N) = 0.50*ac(N)

c-------------------------------
c Print the Fourier coefficients
c-------------------------------

      write (6,*)
      write (6,*) " Fourier coefficients:"
      write (6,*)
      write (6,*) " p, a(p), b(p), |c(p|, a_cs(p), a_sn(p)"
      write (6,*) " --------------------------------------"
      write (6,*)

      Do p=0,NFC
        ps = Dsqrt(af(p)**2+bf(p)**2)
        write (6,101) p,af(p),bf(p),ps,ac(p),as(p)
      End Do

      write (6,*)
      write (6,*)

c--------------------------
c confirm the interpolation
c--------------------------

      write (6,*)
      write (6,*) " Verified values:"
      write (6,*)
      write (6,*) "  i, x, f, F, Fc, Fs"
      write (6,*) " -------------------"

      Do i=1,N1

        xh = x(i)-a

c---
c complete series
c M terms are retained
c---

        z(i) = 0.5*af(0)

        Do p=1,M
         arg  = p*k*xh
         cs   = Dcos(arg)
         sn   = Dsin(arg)
         z(i) = z(i)+af(p)*cs+bf(p)*sn
        End Do

c---
c cosine and sine series
c N terms are retained
c---

        w(i) = 0.5*ac(0)
        v(i) = 0.0

        Do p=1,N
         arg = 0.50D0*p*k*xh
         cs  = Dcos(arg)
         sn  = Dsin(arg)
         w(i)= w(i)+ac(p)*cs
         v(i)= v(i)+as(p)*sn
        End Do

        write (6,101) i,x(i),f(i),z(i),w(i),v(i)

      End Do

c----------
c  Plotting
c----------

      xx  = a-L

      write (6,*) 
      write (6,*) " Continue for plotting ?"
      write (6,*) 
      write (6,*) " Enter 0 for no"
      write (6,*) "       1 for yes"
      write (6,*) " ---------------"

      read (5,*) Icon

      If(Icon.eq.0) Go to 98

      np  = 200
      np1 = np+1
      dx = 3.0*L/np

      write (7,101) np1

      Do 1 i=1,np1

       xh  = xx-a       

       zcm = 0.5*af(0)    ! complete series

       Do p=1,M
         arg = p*k*xh
         cs  = Dcos(arg)
         sn  = Dsin(arg)
         zcm = zcm+af(p)*cs+bf(p)*sn
       End Do

       zcs = 0.5*ac(0)    ! cosine and sine series
       zsn = 0.

       Do p=1,N      
         arg = 0.50D0*p*k*xh
         cs  = Dcos(arg)
         sn  = Dsin(arg)
         zcs = zcs+ac(p)*cs
         zsn = zsn+as(p)*sn
       End Do

c---
c exact value for comparison
c---

       if(menu.eq.1) ww = 1.0
       if(menu.eq.2) ww = xx**2
       if(menu.eq.3) ww = exp(xx)
       if(menu.eq.4) ww = exp(cos(pi2*xx))
       if(menu.eq.5) ww = exp(sin( pi*xx))

       write (6,101) i,xx,ww,zcm,zcs,zsn
       write (7,101) i,xx,ww,zcm,zcs,zsn

       xx = xx+dx

   1  Continue

c----------------
c back to the menu
c----------------

      Go to 98

c-----
c done
c-----

 99   Continue

      write (7,101) Null
      close (7)

 100  Format (1x,f10.5,1x,f10.5,1x,f10.5)
 101  Format (1x,I3,10(1x,f10.5))

      Stop
      End
