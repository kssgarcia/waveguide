      program quad

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c This program accompanies the book:
c
c            C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c------------------------------
c roots of a quadratic equation
c with real coefficients a, b, c
c------------------------------

      Implicit Double Precision (a-h,o-z)

      write (6,*)
      write (6,*) " Roots of: a*x^2 + b*x + c = 0"

 98   Continue

      write (6,*)
      write (6,*) "Please enter the coefficients a, b, c"
      write (6,*) "99 for any one to quit"
      write (6,*) "----------------------"
      read  (5,*) a,b,c 

      if(a.eq.99) Go to 99
      if(b.eq.99) Go to 99
      if(c.eq.99) Go to 99

c------
c traps
c------

      if(abs(a).lt.0.0000001) then

         if(abs(b).gt.0.0000001) then
           write (6,*) 
           write (6,*) " The coefficient a is too small"
           write (6,*) " and the equation looks linear"
           root = -c/b
           write (6,*) " The root is: ", root
         else
           write (6,*) " No root for the constant quadratic"
         end if

         Go to 98

      end if

c-------------------------------------
c proceed with the analytical solution
c-------------------------------------

      d  = b*b-4.0D0*a*c   ! discriminant
      a2 = 2.0D0*a

      if(d.ge.0) then

         srd = Dsqrt(d)
         x1=(-b+srd)/a2
         x2=(-b-srd)/a2

         write (6,*)
         write (6,*) 'roots: ',x1,x2

      else

         d = -d
         real = -b/a2
         rmag = Dsqrt(d)/a2

         write (6,*)
         write (6,*) 'Roots are complex conjugate'
         write (6,*)
         write (6,*) 'Real part:      ',real
         write (6,*) 'Imaginary part: ',rmag

      end if

      Go to 98

c-----
c Done
c-----

  99  Continue

      stop
      end 
