      program inv_u_dr

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c     Oxford University Press
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension a(100,100),b(100,100),v(100,100)

c----------------
c read the matrix
c----------------

      open (2,file="matrix_u.dat")

      read (2,*) N
      Do i=1,N
       read (2,*) (A(i,j),j=i,N)
      End Do

      close (2)

c-----------------
c print the matrix
c-----------------

      write (6,*)
      write (6,*) " Original matrix"
      write (6,*)

      Do i=1,N
       write (6,100) (a(i,j),j=1,N)
      End Do

c--------------------
c compute the inverse
c--------------------

      call inv_u (n,a,b)

c-----------------
c print the inverse
c-----------------

      write (6,*)
      write (6,*) " Alleged inverse"
      write (6,*)

      Do i=1,n
       write (6,100) (b(i,j),j=1,N)
      End Do

c-------
c verify
c-------

      Do i=1,N
        Do j=1,N
         v(i,j) = 0.0D0
         Do k=1,N
          v(i,j) = v(i,j) + a(i,k)*b(k,j)
         End Do
        End Do
      End Do

      write (6,*)
      write (6,*) " Product of A and A(inv) is"
      write (6,*)

      Do i=1,n
       write (6,100) (v(i,j),j=1,n)
      End Do

c-----
c Done
c-----

 100  Format (20(1x,f10.5))

      Stop
      End
