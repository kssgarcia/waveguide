      program bits

c=============================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=============================================

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c---------------------------------------
c Greatest integer that can be described
c with n bits
c
c A bit is a BInary digiT
c A byte consists of 8 bits
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      Integer p,q

c------
c start
c------

      write (6,*)
      write (6,*) " Will compute the greatest integer "
      write (6,*) " that can be described with n bits "
      write (6,*) " --------------------------------- "

 98   Continue        ! return to repeat

      write (6,*) 
      write (6,*) " Enter the number of bits (less than 32)"
      write (6,*) " 0 to quit "
      write (6,*) " ---------------"
      read  (5,*) n

      If(n.eq.0) Go to 99

      write (6,*)
      write (6,*) " bits     increment       largest integer"
      write (6,*)

      q = 0

      Do i=0,n-1
        p = 2**i
        q = q+p
        write (6,100) i+1,p,q
      End Do

      Go to 98

c-----
c done
c-----

  99  continue

 100  Format (1x,i5,2(1x,i15))

      stop
      end
