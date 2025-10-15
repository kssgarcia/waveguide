      Program golden_ratio

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c--------------------------------------
c Computes the golden ratio in terms of 
c the Fibonacci series
c--------------------------------------

      Implicit Double precision (a-h,o-z)

      Dimension a(512)

      Parameter(ncount=10)

c---------
c start up
c---------

      a(1) = 1
      a(2) = 1
      k = 3

      Icount = 1
  98  Continue

      a(k) = a(k-1)+a(k-2)
      ratio = a(k)/a(k-1)

      write (6,100) k,ratio

c---------------
c reset counters
c---------------

      k = k+1

      icount = icount + 1

      If(icount.eq.ncount) then
        icount = 1
        write (6,*) 
        write (6,*) " More numbers ? "
        write (6,*) " 0 for NO, 1 for YES"
        write (6,*) " -------------------" 
        write (6,*) 
        read  (5,*) more
        If(more.eq.0) Go to 99
      End If

      Go to 98

c-----
c Done
c-----

 99   Continue

 100  Format (1x,i4,2x,f10.8)

      Stop
      End
