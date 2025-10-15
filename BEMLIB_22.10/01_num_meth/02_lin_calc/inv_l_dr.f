      program inv_l_dr

c-----------------------------------------
c Copyright 1999 by C. Pozrikidis
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c       Oxford University Press
c------------------------------------------------

c-------------------------------------
c Inverse of a lower triangular matrix
c-------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension a(128,128),b(128,128),v(128,128)

c----------------
c read the matrix
c----------------

      open (2,file="matrix_l.dat")

       read (2,*) N
       Do i=1,N
        read (2,*) (a(i,j),j=1,i)
       End Do

      close (2)

c-----------------
c print the matrix
c-----------------

      write (6,*)
      write (6,*) " Original matrix:" 
      write (6,*)

      Do i=1,N
       write (6,100) (a(i,j),j=1,i)
      End Do

c--------------------
c compute the inverse
c--------------------

      call inv_l (n,a,b)

c-----------------
c print the inverse
c-----------------

      write (6,*)
      write (6,*) " Alleged inverse:" 
      write (6,*)

      Do i=1,N
       write (6,100) (b(i,j),j=1,i)
      End Do

c---
c Verify
c---

      Do i = 1,N
        Do j = 1,N
         v(i,j) = 0.0D0
         Do k=1,N
          v(i,j) = v(i,j) + a(i,k)*b(k,j)
         End Do
        End Do
      End Do

      write (6,*)
      write (6,*) " The product of A and A(inv) is:"
      write (6,*)

      Do i=1,n
       write (6,100) (v(i,j),j=1,n)
      End Do

c-----
c Done
c-----

 100  Format (20(1x,f8.5))

      stop
      end
