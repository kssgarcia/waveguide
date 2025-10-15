      program gauss_log_dr 

c==========================================
c          FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press
c------------------------------------------------

c----------------------------------------
c Compute the integral of a function with a logarithmic
c singularity using the Gauss-log quadrature
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension z(20),w(20),y(20)

      common/pii/pi

c----------
c constants
c----------

      pi = 3.14159 265358 D0

c------------
c preferences
c------------
 
      open (2,file="gauss_log.out")

      write (6,*) 
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 1 for q(x) = 1.0'
      write (6,*) ' 2 for q(x) = example in text'
      write (6,*)
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c----
c Number of points
c---

  97  Continue

      write (6,*)
      write (6,*) ' Enter the number of base points'
      write (6,*) ' Select from 1, 2, 3, 4, 5'
      write (6,*)
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) NQ

      if(NQ.eq.0) Go to 99

c----------------------
c abscissae and weights
c----------------------

      call gauss_log (NQ,z,w)

c----------------------
c function evaluations
c----------------------

      Do i=1,NQ
        call gauss_log_fun (z(i),y(i),menu)
      End Do

c----------------------
c computation of the integral
c----------------------

      sum = 0.0D0

      Do i=1,NQ
        sum = sum+y(i)*w(i)
      End Do

c--------
c display
c--------

      write (6,102) NQ,sum
      write (6,*) 

      Go to 97 ! repeat

c-----
c Done
c-----

  99  Continue

  100 Format (1x,I3,1x,I5,3(1x,f15.10))
  101 Format (1x,I3,3(1x,f15.10))
  102 Format (" Number of base points:",I3," Integral:",f15.10)

      Stop
      End
