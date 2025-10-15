      program lobatto_dr

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c-----------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c       Oxford University Press
c-----------------------------------------------------

c------------------------------------------------
c compute the integral of a nonsingular function 
c using the Lobatto quadrature
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  z(20), w(20)
      Dimension  f(20)

      common/pii/pi

c----------
c constants
c----------

      pi = 3.14159 265358 D0

c------------
c preferences
c------------
 
      open (2,file="lobatto.out")

      write (6,*) 
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 1 for f(x) = 1.0'
      write (6,*) ' 2 for f(x) = exponential'
      write (6,*) ' 3 for f(x) = x^4'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c-----------------
c Number of points
c-----------------

  97  Continue

      write (6,*) ' Please enter the number of base points (NQ)'
      write (6,*) ' Choose from 2,3,4,5,6,7, 8'
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) NQ

      if(NQ.eq.0) Go to 99

c----------------------
c abscissae and weights
c----------------------

      call lobatto (NQ,z,w)

c---------------------
c function evaluations
c---------------------

c     write (6,101) k+1

      Do i=1,NQ
        call lobatto_dr_fun (z(i),f(i),menu)
c       write (6,101) i,z(i),f(i),w(i)
      End Do

c----------------------------
c computate the integral
c----------------------------

      sum = 0.0D0

      Do i=1,NQ
        sum = sum+f(i)*w(i)
      End Do

c--------
c display
c--------

      write (6,*) 
      write (6,102) NQ,sum
      write (6,*) 

      Go to 97 ! repeat

c-----
c done
c-----

  99  Continue

  100 Format (1x,I3,1x,I5,3(1x,f15.10))
  101 Format (1x,I3,3(1x,f20.15))
  102 Format (" Base points (NQ): ",I3," Integral:",f15.10)

      stop
      end
