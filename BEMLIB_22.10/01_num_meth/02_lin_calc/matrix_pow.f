      program matrix_pow

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c------------------------------------------------

c-------------------
c Powers of a matrix
c-------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),b(128,128),c(128,128)

c----------------
c read the matrix
c----------------

      open (2,file="matrix.dat")

      read (2,*) N
      Do i=1,N
       read (2,*) (a(i,j),j=1,N)
      End Do

      close (2)

c-----------------
c print the matrix
c-----------------

      write (6,*)
      write (6,*) " Matrix A:" 
      write (6,*)

      Do i=1,N
       write (6,100) (a(i,j),j=1,N)
      End Do

c--------------------------------------
c initialize the current power matrix b
c--------------------------------------

      Do i=1,N
       Do j=1,N
        b(i,j) = a(i,j)
       End Do
      End Do

      Iexp = 2

c-----------------------
c compute the next power
c-----------------------

  98  Continue

      Do i=1,N
       Do j=1,N
       c(i,j) = 0.0
        Do k=1,N
         c(i,j) = c(i,j)+a(i,k)*b(k,j)
        End Do
       End Do
      End Do

c------
c print
c------

      write (6,*)
      write (6,*) Iexp,"-power of the matrix A: " 
      write (6,*)

      Do i=1,N
       write (6,100) (c(i,j),j=1,N)
      End Do

      write (6,*)
      write (6,*) " One more power? "
      write (6,*) " 1 for yes, 0 for no"
      read  (5,*) more

      if(more.eq.0) Go to 99

      Iexp = Iexp+1

c----
c replace the old with the new power
c----

      Do i=1,N
       Do j=1,N
          b(i,j) = c(i,j)
       End Do
      End Do

      Go to 98

c-----
c Done
c-----

  99  Continue

 100  Format (20(1x,f12.8))

      Stop
      End
