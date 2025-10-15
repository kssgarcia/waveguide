      program integral_1d

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c====================================================
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c====================================================

c-----------------------------------------------------
c Computation of a one-dimensional
c integral of a non-singular function
c by one of the following methods:
c
c (a) Trapezoidal rule
c (b) Simpson's 1/3 rule
c (c) Gauss-Legendre quadrature
c-----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension x(90000),y(90000)
      Dimension Z(90000),W(90000)

      Dimension ZZ(20),WW(20)

      common/ccc/rks,cnt
      common/vari/m
      common/varr/c
      common/pii/pi
      common/option72/alpha72,p72,q72
      common/option73/alpha73,pmq73

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi
      pi2 = 2.0D0*pi
      pi4 = 4.0D0*pi

      null= 0

      tpf = 2.5D0

c--------
c prepare
c--------
 
      open (2,file="integral_1d.out")

c------------------
c menu of functions
c------------------

      write (6,*) 
      write (6,*) ' Menu of functions'
      write (6,*) 
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) '  1 for f(x) = -2 x^2 (lnx - 1) exp(-x^2)'
      write (6,*) '  2 for f(x) = 3.0'
      write (6,*) '  3 for f(x) = ln(1-0.99 cos^2(x*pi/2))'
      write (6,*) '  4 for f(x) = exp(x)                    '
      write (6,*) '  5 for f(x) = x**3/2                   '
      write (6,*) '  6 for f(x) = exp(-(1+t)**2/4)         '
      write (6,*) '  7 for f(x) = (1-exp(-t))/t^(3/2)      '
      write (6,*) ' 20 for the complete second kind elliptic integral'
      write (6,*) ' 21 for f(x) = ln|c-sin(x)|              '
      write (6,*) ' 22 for f(x) = 2*(1-exp(-t^2/4))/t^2      '
      write (6,*) ' 30 for f(x) = (1-cos(mt))/sin(t/2)**2    '
      write (6,*) ' 31 for f(x) = 1/(m**2+sin(t/2)**2)    '
      write (6,*) ' 32 for f(x) = t**2/sqrt(2-cos(t)**2)    '
      write (6,*) ' 41 for f(x) = ( cos(2mt)-cos(2(m-1)t) )/sint  '
      write (6,*) ' 42 for f(x) = ( cos(2mt)-cos(2(m-1)t) ) '
      write (6,*) '                /(sint*sqrt(1+2sin^2))  '
      write (6,*) ' 45 for f(x) = -(1/pi) ( 1-cos(2mt) ) '
      write (6,*) '               /sqrt(cos(2x)**2-8*cos(2x)+7)  '
      write (6,*) ' 60 for f(x) = cos(cx) * (1-x**2)   '
      write (6,*) ' 61 for f(x) = (sinx/x-cosx)/x**2'
      write (6,*) ' 70 for f(x) = exp(-0.5*sinx)'
      write (6,*) ' 71 for f(x) = exp(-2.0*sinx)'
      write (6,*) ' 72 for f(x) = sin(x/2)^alpha*sin(px)*sin(qx)'
      write (6,*) ' 73 for f(x) = sin(x/2)^alpha*cos((p-q)x)'
      write (6,*) 
      write (6,*) ' 80 for f(x) = complicated1'
      write (6,*) 
      write (6,*) ' 0 to quit'
      write (6,*) ' ----------'
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c------------------
c choose the method
c------------------

      write (6,*)
      write (6,*) ' Choose the integration method'
      write (6,*)
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 1 for trapezoidal rule'
      write (6,*) ' 2 for Simpson 1/3 rule'
      write (6,*) ' 3 for Gauss--Legendre quadrature'
      write (6,*) 
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) method

      if(method.eq.0) Go to 99

c-----------------
c return to repeat
c-----------------

      Ipass = 0
   97 Continue
      Ipass = Ipass + 1

c--------
c prepare
c--------

      if(menu.eq.20) then

        a = 0.0
        b = pih
        write (6,*)
        write (6,*) ' Please enter the value of k'
        write (6,*) ' ---------------------------'
        read  (5,*) rk
        rks = rk*rk

      elseif(menu.eq.21) then

        a = -pih
        b =  pih

        write (6,*)
        write (6,*) ' Please enter the constant: -1<c<1'
        write (6,*) ' 99 to quit '
        write (6,*) ' -----------'
        read  (5,*) cnt
        if(abs(cnt-99.0).lt.0.000001) Go to 99

      elseif(menu.eq.30) then
        write (6,*)
        write (6,*) ' Please enter m'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) m
        if(m.eq.0) Go to 99
        a = 0.0001D0
        b = 0.9999*pi2

      elseif(menu.eq.31) then

        write (6,*)
        write (6,*) ' Please enter m'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) m
        if(m.eq.0) Go to 99
        a = 0.0001D0
        b = 0.9999*pi2

      elseif(menu.eq.32) then

        a = -1.0D0
        b =  1.0D0

      elseif(menu.eq.41.or.menu.eq.42) then

        write (6,*)
        write (6,*) ' Please enter m'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) m
        if(m.eq.0) Go to 99
        a = 0.00001D0
        b = 0.9999*pi
        if(menu.eq.45) b = 0.5D0*pi;

      elseif(menu.eq.45) then

        write (6,*)
        write (6,*) ' Please enter m'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) m
        if(m.eq.0) Go to 99
        a = 0.00001D0
        b = 0.50000D0*pi;

      elseif(menu.eq.60) then

        write (6,*)
        write (6,*) ' Please enter the constant c'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) c
        if(c.eq.0) Go to 99
        a = -1.0D0
        b =  1.0D0

      elseif(menu.eq.61) then

        a = 0.00001D0
        b = 100.0D0

      elseif(menu.eq.70) then

        a = 0.00001D0
        b = 2.0*pi

      elseif(menu.eq.71) then
        a = 0.00001D0
        b = 2.0*pi

      elseif(menu.eq.72) then

        a = 0.00001D0
        b = pi
        alpha72 = 1.0D0
        write (6,*) ' Please enter p and q'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) p72,q72
        if(p72.eq.0) Go to 99
        if(q72.eq.0) Go to 99

      elseif(menu.eq.73) then

        a = 0.00001D0
        b = pi
        write (6,*) ' Please enter alpha and p-q'
        write (6,*) ' 0 to quit'
        write (6,*) ' ---------'
        read  (5,*) alpha73,pmq73

      elseif(menu.eq.80) then
        a = 0.0000D0
        b = pi/3.0D0

      else
        write (6,*)
        write (6,*) ' Enter the lower and upper limits a and b'
        write (6,*) ' ----------------------------------------'
        read  (5,*) a,b

      end if

c--------------------------
c exact value, if available
c--------------------------

      Iexact = 0

      if(menu.eq.4) then
         exact = exp(b)-exp(a)
         Iexact = 1

      elseif(menu.eq.5) then
         exact = (b**tpf-a**tpf)/tpf
         Iexact = 1

      elseif(menu.eq.20) then

         if(rk.eq.0.8) then
           exact = 1.276349943
         elseif(rk.eq.0.9) then
           exact = 1.171697053
         else
           write (6,*)
           write (6,*) ' Please enter exact value (if known)'
           write (6,*) ' -----------------------------------'
           read  (5,*) exact
         end if

         Iexact = 1

      elseif(menu.eq.21) then

         exact = -pi*log(2.0)
         Iexact = 1

      elseif(menu.eq.22) then

         exact = sqrt(pi)
         Iexact = 1

      elseif(menu.eq.30) then
         exact = 0.5*m
         Iexact = 1

      elseif(menu.eq.31) then
         exact = 2.0*pi/m/sqrt(m**2+1.0)
         Iexact = 1

      elseif(menu.eq.32) then
         exact = 0.5*pi-1.0D0
         Iexact = 1

      elseif(menu.eq.41) then
         exact = -1.0D0/pi/(2.0D0*m-1.0D0)
         Iexact = 1

      elseif(menu.eq.60) then
         exact = 4.0/c**2 *(sin(c)/c-cos(c))
         Iexact = 1

      elseif(menu.eq.61) then
         exact = pi/4.0D0;
         Iexact = 1

      elseif(menu.eq.80) then
         exact = 3.0/160.0;
         Iexact = 1

      end if

c---------------------
c compute the integral
c---------------------

c---
c return to repeat
c---

      Jpass = 0
   98 Continue
      Jpass = Jpass + 1

c---
c number of base points
c---

      if(method.eq.1) then

         write (6,*)
         write (6,*) ' Enter the number of intervals N'
         write (6,*)
         write (6,*) ' 0 to quit '
         write (6,*) ' ----------'
         read  (5,*) N

         If(N.eq.0) Go to 99

         N1 = N+1

      elseif(method.eq.2) then

         write (6,*)
         write (6,*) 'Enter the number of intervals N (even)'
         write (6,*)
         write (6,*) '0 to quit'
         write (6,*) '---------'
         read  (5,*) N

         If(N.eq.0) Go to 99

         N1 = N+1

      else

         write (6,*)
         write (6,*) 'Enter the number of base points'
         write (6,*) 'choose from: 1,2,3,4,5,6,8,12,20'
         write (6,*)
         write (6,*) '0 to quit '
         write (6,*) '----------'
         read  (5,*) NGL

         if(NGL.eq.0) Go to 99

         call gauss_leg (NGL,ZZ,WW)

         Do i=1,NGL
           Z(i) = ZZ(i)
           W(i) = WW(i)
         End Do

      End If

c--------------------
c integration weights
c--------------------

      if(method.eq.1) then            ! Trapezoidal rule

         w(1) = 0.5D0
         Do i=2,N
           w(i) = 1.0D0
         End Do
         w(N1)= 0.5D0

      elseif(method.eq.2) then       ! Simpson's rule

         w(1) = 1.0D0
         Do i=2,N,2
           w(i)   = 4.0D0
           w(i+1) = 2.0D0
         End Do
         w(N1)= 1.0D0

      End If

c---------------------
c integrand evaluation
c and computation of the integral
c---------------------

c---
      if(method.le.2) then
c---

        h = (b-a)/N

        Do i=1,N1
          x(i) = a+h*(i-1.0D0)
          call integral_1d_fun(x(i),y(i),menu)
        End Do

       sum = 0.0D0
       Do i=1,N1
        sum = sum+y(i)*w(i)
       End Do

       if(method.eq.1) sum = sum*h
       if(method.eq.2) sum = sum*h/3.0D0

       nbp = N1

c---
      elseif(method.eq.3) then
c---

        xm = 0.5D0*(b+a)
        xd = 0.5D0*(b-a) 

       Do i=1,NGL
          x(i) = xm+xd*z(i)
          call integral_1d_fun(x(i),y(i),menu)
       End Do

       sum = 0.0D0
       Do i=1,NGL
         sum = sum+y(i)*w(i)
       End Do

       sum = sum*xd

       nbp = NGL

c---
      End If
c---

c----------------------------

       if(menu.eq.30) sum =  sum/pi4
       if(menu.eq.41) sum =  sum/pi4
       if(menu.eq.42) sum =  sum/pi4
       if(menu.eq.45) sum = -sum/pi
       if(menu.eq.72) sum = 2**(alpha72+1.0)*sum/pi
       if(menu.eq.73) sum = 2**(alpha73+0.0)*sum/pi

       write (6,102) nbp,sum

c------------------
c compute the error
c-----------------

      if(Iexact.eq.1) then

        error = sum-exact
        ord   = log(abs(error))
        absi  = log(h)

        write (2,101) Jpass,h,sum,error,absi,ord
        write (6,*) 
        write (6,108) 
        write (6,*) 
        write (6,191) Jpass,h,sum,exact,error,absi,ord
        write (6,*) 

      end If
 
c-------------------------------
c return for another computation
c-------------------------------

      Go to 97
c     Go to 98

c--------
c wrap up
c--------

  99  Continue

      write (2,100) null
      close (2)

c-----
c Done
c-----

  100 Format (1x,I3,1x,I5,3(1x,f15.10))
  101 Format (1x,I3,8(1x,f13.7))
  191 Format (1x,I3,9(1x,f13.7))
  102 Format (1x,"Number of points :",I6," Integral",f20.15)
  108 Format ("trial, step, integral, exact,error"
     +       ,", log(step), log(error)")

      Stop
      End
