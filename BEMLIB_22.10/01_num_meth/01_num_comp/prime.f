      program prime

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c       C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

c----------------------------------------
c Checks whether a given integer is prime
c----------------------------------------

 98   Continue

      write (6,*)
      write (6,*) " Please enter the integer to be tested"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read (5,*) n

      If(n.eq.0) Go to 99

      If(n.le.0) then
        write (6,*)
        write (6,*) " The integer must be positive; try again"
        write (6,*)
        Go to 98
      End If

c--------
c testing
c--------

      Do m=2,n-1

        l = n/m    ! testing for the remainder
        k = l*m

        If(k.eq.n) then
          write (6,*) 
          write (6,*) n,' is not a prime number'
          write (6,*) 
          write (6,*) ' Its highest divisor is ',l
          write (6,*)
          Go to 98
        End If 

      End Do

      write (6,*) n,' is a prime number'

      Go to 98   ! Return to repeat

c-----
c Done
c-----

 99   Continue

      Stop
      End
