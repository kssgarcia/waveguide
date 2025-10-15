      program matrix_mult

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c              C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c------------------------------------------------

c------------------------------------------------
c  Multiply two square matrices
c
c  SYMBOLS:
c  -------
c
c  n...dimension of the matrices
c  a...first matrix
c  b...second matrix
c
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(50,50),b(50,50),c(50,50)

c----------------
c read the matrix
c and the vector
c---------------

      write (6,*)
      write (6,*) " Reading the matrices from file: matrix_ab.dat "
      write (6,*) " ----------------------------------------------"
      write(6,*)

      open (unit=8,file="matrix_ab.dat",status="unknown")

      read (8,*) n
      read (8,*)

      Do i=1,n
       read  (8, *)  (a(i,j),j=1,n)
       write (6,100) (a(i,j),j=1,n)
      End Do

      read (8,*)
      write(6,*)

      Do i=1,n
       read  (8,*)  (b(i,j),j=1,n)
       write (6,100) (b(i,j),j=1,n)
      End Do

      close (8)
 
c---------
c multiply
c---------

      Do i=1,n
       Do j=1,n
        c(i,j) = 0.0D0
        Do l=1,n
         c(i,j) = c(i,j) + a(i,l)*b(j,j)
        End Do
       End Do
      End Do

c--------
c Display
c--------


      write (6,*)
      write (6,*) " Product:"
      write (6,*) " --------"

      Do i=1,n
       write (6,100) (c(i,j),j=1,n)
      End Do

c-----
c Done
c-----

 100  Format (5(2x,f10.5))
 101  Format (20(2x,f10.5))

      Stop
      End
