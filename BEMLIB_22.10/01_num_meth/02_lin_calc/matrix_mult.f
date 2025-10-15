      program matrix_mult

c-----------------------------------------
c Copyright 1999 by C. Pozrikidis
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c------------------------------------------------

c---------------------------------
c Multiply two matrices: c = a . b
c---------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension A(128,128),B(128,128),C(128,128)

c------------------
c read the matrices
c------------------

      open (2,file="matrix2.dat")

      read (2,*) M,N
      Do i=1,M
       read (2,*) (A(i,j),j=1,N)
      End Do

      read (2,*) N,K
      Do i=1,N
       read (2,*) (B(i,j),j=1,K)
      End Do

      close (2)

c-----------------
c print the matrices
c-----------------

      write (6,*)
      write (6,*) " Matrix A:" 
      write (6,*)

      Do i=1,M
       write (6,100) (A(i,j),j=1,N)
      End Do

      write (6,*)
      write (6,*) " Matrix B:" 
      write (6,*)

      Do i=1,N
       write (6,100) (B(i,j),j=1,K)
      End Do

c----------------------------
c compute the product c = a*b
c----------------------------

      Do i=1,M
       Do j=1,K
        C(i,j) = 0.0D0
        Do l=1,N
         C(i,j) = C(i,j) + A(i,l)*B(l,j)
        End Do
       End Do
      End Do

c------------------
c print the product
c------------------

      write (6,*)
      write (6,*) " Matrix C:" 
      write (6,*)

      Do i=1,M
       write (6,100) (C(i,j),j=1,K)
      End Do


c-----
c Done
c-----

  99  Continue

 100  Format (20(1x,f12.8))

      Stop
      End
