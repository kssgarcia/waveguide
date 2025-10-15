      program euclid

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c          Oxford University Press
c------------------------------------------------

c---------------------------------
c This program finds the greatest common divisor
c of two positive integers using Euclid's algorithm
c---------------------------------

      write (6,*)
      write (6,*) " Will compute the Gratest Common Divisor"
      write (6,*) "      of two positive integers"
      write (6,*)

 98   write (6,*)
      write (6,*) " Please enter the two integers"
      write (6,*) "      0 for either one to quit"
      write (6,*) " -----------------------------"
      read  (5,*) n,m

      If(n.eq.0) Go to 99
      If(m.eq.0) Go to 99

  1   Continue

      If(n.eq.m) then
        k=n
        Go to 2
      End If

      If(n.gt.m) then
        n1 = m
        m  = n
        n  = n1
      End If

      k = m-n
      m = n
      n = k

      Go to 1

   2  Continue

      write (6,*)
      write (6,100) k

      Go to 98   ! repeat

c-----
c Done
c-----

  99  Continue

      write (6,*) " Thank you for running me; run me again"

 100  Format (" The Greatest Common Divisor is: ",i7)

      Stop
      End
