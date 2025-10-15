      program gauss_trgl_dr 

c============================================
c          FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c============================================

c------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press
c------------------------------------------------

c-----------------------------------------------
c Computes the integral of a regular function 
c over the standard triangle in the xi-eta plane
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension xi(20),eta(20),w(20),f(20)

      common/pii/pi

c----------
c constants
c----------

      pi = 3.14159 265358 D0

c------------
c preferences
c------------
 
      open (2,file="gauss_trgl.out")

      write (6,*) 
      write (6,*) ' Enter:'
      write (6,*) 
      write (6,*) ' 1 for f(xi, eta) = 1.0'
      write (6,*) ' 2 for f(xi, eta) = xi'
      write (6,*) ' 3 for f(xi, eta) = eta'
      write (6,*) ' 4 for f(xi, eta) = xi^2'
      write (6,*) ' 5 for f(xi, eta) = eta^2'
      write (6,*) ' 6 for f(xi, eta) = xi*eta'
      write (6,*) ' 7 for f(xi, eta) = xi^3'
      write (6,*) ' 8 for f(xi, eta) = xi^2*eta'
      write (6,*)
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) menu

      if(menu.eq.0) Go to 99

c----
c number of points
c---

  97  Continue

      write (6,*)
      write (6,*) ' Enter the number of base points'
      write (6,*) ' Select from 1,3,4,6,7,9,12,13'
      write (6,*)
      write (6,*) ' 0 to quit'
      write (6,*) ' ---------'
      read  (5,*) mint

      if(mint.eq.0) Go to 99

c----------------------
c Abscissae and weights
c----------------------

      call gauss_trgl (mint,xi,eta,w)

c---------------------
c function evaluations
c---------------------

      Do i=1,mint
        call gauss_trgl_fun (xi(i),eta(i),f(i),menu)
      End Do

c----------------------
c computate the integral
c----------------------

      sum = 0.0D0

      Do i=1,mint
        sum = sum+f(i)*w(i)
      End Do

      sum = 0.50D0*sum

c--------
c display
c--------

      write (6,102) mint,sum
      write (6,*) 

      Go to 97 ! repeat

c-----
c Done
c-----

  99  Continue

  100 Format (1x,I3,1x,I5,3(1x,f15.10))
  101 Format (1x,I3,3(1x,f15.10))
  102 Format (" Number of quad points (NQ):",I3," Integral:",f15.10)

      Stop
      End
