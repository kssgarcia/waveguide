      program gauss_lag_dr

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c----------------------------------------------------
c This program accompanies the book:
c
c                 C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c            Oxford University Press
c----------------------------------------------------

c----------------------------------------------------
c Compute the integral of an exponentially decaying
c function using the Gauss--Laguerre quadrature
c----------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension z(2000),w(2000),y(2000)

      common/pii/pi

c----------
c constants
c----------

      pi = 3.14159 265358 D0

c------------
c preferences
c------------
 
      write (6,*)
      write (6,*) '        MENU'
      write (6,*)
      write (6,*) ' Enter 1 for q(x) = 1.0 '
      write (6,*) '       2 for q(x) = x '
      write (6,*) '       3 for q(x) = x^2 '
      write (6,*) '       4 for q(x) = x^3 '
      write (6,*) '       5 for q(x) = x^4 '
      write (6,*) '       6 for q(x) = x^5 '
      write (6,*) '       7 for q(x) = x^6 '
      write (6,*) '       8 for q(x) = x^7 '
      write (6,*) '       9 for q(x) = x^8 '
      write (6,*) '      10 for q(x) = x^9 '
      write (6,*) '      11 for q(x) = x^10 '
      write (6,*) '      20 for example in the text'
      write (6,*) '      25 another example '
      write (6,*)
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c--------
c prepare
c--------

      open (2,file="integral_lag.out")

c-----------------
c number of points
c-----------------

  97  Continue

      write (6,*) 
      write (6,*) ' Enter NQ'
      write (6,*) ' choose from 1, 2, 3, 4, 5'
      write (6,*)
      write (6,*) ' 0 to quit'
      write (6,*) ' ----------'
      read  (5,*) NQ

      if(NQ.eq.0) Go to 99

c----------------------
c abscissae and weights
c----------------------

      call gauss_lag (NQ,z,w)

c----------------------
c evaluate the function
c----------------------

      Do i=1,NQ
        call gauss_lag_fun (z(i),y(i),menu)
      End Do

c---------------------
c compute the integral
c---------------------

      sum = 0.0D0

      Do i=1,NQ
        sum = sum + y(i)*w(i)
      End Do

      write (6,102) NQ,sum
      write (6,*) 

c----
c repeat
c----

      Go to 97

c-----
c Done
c-----

  99  Continue

  100 Format (1x,I3,1x,I5,3(1x,f15.10))
  101 Format (1x,I3,3(1x,f15.10))
  102 Format (1x,"Number of quad points (NQ):",I3," Integral",f15.6)

      Stop
      End
