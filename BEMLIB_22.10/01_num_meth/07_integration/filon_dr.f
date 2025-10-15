      program filon_dr

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
c Numerical Computation in Science and Engineering
c         Oxford University Press
c------------------------------------------------

c----------------------------------------
c Driver for computing the integral of
c an oscillatory integrand
c using Filon's method
c
c (Pozrikidis 1998), pp. 369-370
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision k

c------
c input
c------

      write (6,*) 
      write (6,*) '   MENU OF FUNCTIONS'
      write (6,*) 
      write (6,*) 'Enter:'
      write (6,*) 
      write (6,*) '1 for q(x) = 1'
      write (6,*) '2 for q(x) = x'
      write (6,*) '3 for q(x) = x^2'
      write (6,*) '4 for q(x) = x^3'
      write (6,*) '5 for q(x) = exp(-x)'
      write (6,*) 
      write (6,*) '0 to quit'
      write (6,*) '---------'
      read(5,*) menu

      if(menu.eq.0) Go to 99

      write (6,*) 
      write (6,*) 'Please enter the lower limit a'
      write (6,*) '             the upper limit b'
      write (6,*) '             the wave number k'
      write (6,*) '------------------------------'
      read  (5,*) a,b,k

 98   Continue

      write (6,*)
      write (6,*) 'Please enter the number of intervals n'
      write (6,*) 'Must be even and less than 902'
      write (6,*)
      write (6,*) '0 to quit'
      write (6,*) '---------'
      read (5,*) n

      if(n.eq.0) Go to 99

c---
c traps
c---

      if(n.gt.900) then
        write (6,*)' unacceptable n > 898; will set n = 898'   
        n = 900
      end if

      if(mod(n,2).eq.1) then
        write (6,*) 
        write (6,*) '<Number of intervals must be even>' 
        Go to 98
      end if

c---
c compute the integrals Ic and Is
c equations (7.6.2)
c---

      call filon 
     +
     +   (menu
     +   ,a,b,k
     +   ,n
     +   ,cos_int,sin_int
     +   )

      write (6,*)
      write (6,110) cos_int
      write (6,100) sin_int
      write (6,*)

      Go to 98      ! return to repeat

c---
c Done
c---

   99 Continue

 110  Format ('cos integral = ', f15.10)
 100  Format ('sin integral = ', f15.10)

      Stop
      End
