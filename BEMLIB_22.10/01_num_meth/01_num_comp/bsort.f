       program ran_bsort

c=========================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c------------------------------------------------

c-------------------------------------
c Generate a set of random numbers
c and then sort them using bubble sort
c-------------------------------------

      Dimension x(3000)

  98  Continue

      write (6,*)
      write (6,*) " How many numbers ?"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) n

      if(n.eq.0) Go to 99

c---------------------------
c Generate the random numbers
c---------------------------

      idum = 1

      write (6,*)
      write (6,*) " Generated numbers:"
      write (6,*) " ------------------"

      Do i=1,n

c       x(i) = r_lcran()       ! on SUN
        x(i) = ran2(idum)
        write (6,100) i,x(i)
      End Do

      write (6,*)

c---------------
c bubble sorting
c---------------

      k = n-1
 1    Continue
      Istop = 1

      Do i=1,k
        i1 = i+1
        if(x(i).gt.x(i1)) then
          save  = x(i)
          x(i)  = x(i1)
          x(i1) = save
          Istop = 0
        End if
      End Do

      k = k-1

      if(Istop.eq.0) Go to 1 

c---------
c printing
c---------

      write (6,*)
      write (6,*) " Sorted numbers:"
      write (6,*) " ---------------"

      Do i=1,N
        write (6,100) i,x(i)
      End Do

      Go to 98      ! repeat

c-----
c Done
c-----

  99  Continue

 100  Format (1x,i4,1x,f12.5,1x,i4)

      Stop
      End
