      program bear_stow

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
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press, 1998
c------------------------------------------------

c---------------------------------
c Driver for
c Bairstow's method for computing
c all roots of a polynomial
c---------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(20),x(20,2)

c----------
c inquiries
c----------

 98   Continue

      write (6,*) 
      write (6,*) " Please enter order of the polynomial n"
      write (6,*) " 0 to quit"
      write (6,*) " ---------" 
      read  (5,*) n

      If(n.eq.0) Go to 99

      n1 = n+1

      write (6,*) 
      write (6,*) " a(1) is the coefficient of x**n"
      write (6,*) 
      write (6,*) " Please enter coefficients starting with a(1)"
      write (6,*) " a total of ",n1," values"
      write (6,*) "-------------------------" 
      read  (5,*) (a(i),i=1,n1)

c------------
c call solver
c------------

      call bairstow (a,n,x,Iflag)

c---
c printing
c---

      write (6,*) 
      write (6,*) " ROOTS:"
      write (6,*) 
      write (6,*) "       Real           Imaginary"
      write (6,*) 

      Do i=1,n
         write (6,100) i,x(i,1),x(i,2)
      End Do

      Go to 98

c-----
c Done
c-----

 99   Continue

  100 Format (1x,i2,1x,f15.10,1x,f15.10)

      Stop
      End
