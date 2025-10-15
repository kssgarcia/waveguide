      program isort

c=======================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c----------------------------------------
c Generate a set of random numbers
c and then index and sort them
c----------------------------------------

c----------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------

      Dimension x(3000)
      Integer   m(3000)

  98  Continue

      write (6,*) 
      write (6,*) " How many random numbers ?"
      write (6,*) " less then 3000 please"
      write (6,*) " 0 to quit"
      write (6,*) " ---------"
      read  (5,*) n

      If(n.eq.0) Go to 99

c---------------------
c Generate the numbers
c---------------------

      write (6,*) 
      write (6,*) " List of random numbers:"
      write (6,*) 

      idum = 1

      Do i=1,n
        x(i) = ran2(idum)
        write (6,100) i,x(i)
      End Do

      write (6,100)

c----------------------------
c indexing
c
c greatest number has index 1
c----------------------------

      write (6,*) 
      write (6,*) " Indexed list:"
      write (6,*) 

      Do i=1,N
         m(i)=1
         Do j=1,N
           if(x(i).lt.x(j)) m(i)=m(i)+1
         End Do
         write (6,100) i,x(i),m(i)
      End Do

      write (6,100)
 
c----------------------------
c Printing in descending order
c----------------------------

       write (6,*) 
       write (6,*) " Sorted list:"
       write (6,*) 

      Do i=1,N
        Do j=1,N
         if(i.eq.m(j)) write (6,100) j,x(j),m(j)
        End do
      End Do

      Go to 98    ! repeat

c-----
c Done
c-----

  99  Continue

 100  Format (1x,i4,1x,f12.5,1x,i4)

      Stop
      End
