      program ldu_dec

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press
c
c 1998
c------------------------------------------------

c------------------------------------------
c LDU decomposition of a symmetric matrix
c     computation of the determinant
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),v(10,10)
      Dimension d(10),u(10,10)
      Double precision l(10,10)
	
c----------------
c Read the matrix
c----------------

      open (unit=3,file='matrix_sym.dat')

      read (3,*) n
      Do i=1,n
        read (3,*) (a(i,j),j=1,n)
      End Do

      close(3)

c---
c printing session
c---

      write (6,*)
      write (6,*) " Original matrix A"
      write (6,*)

      Do i=1,n
        write (6,100) (a(i,j),j=1,n)
      End Do 

c--------------------
c Do the decomposition
c---------------------

      call ldu (n,a,l,d,u)

c---
c printing session
c---

      write (6,*)
      write (6,*) " matrix L"
      write (6,*)

      Do i=1,n
        write (6,100) (l(i,j),j=1,n)
      End Do 

      write (6,*)
      write (6,*) " vector d"
      write (6,*)

      write(6,100) (d(j),j=1,n)

      write (6,*)
      write (6,*) " matrix U"
      write (6,*)

      Do i=1,n
        write (6,100) (u(i,j),j=1,n)
      End Do 
	
c---
c verify
c---

      Do i=1,n
        Do j=1,n
          v(i,j) = 0.0D0
          k = i
          if(j.lt.k) k = j
          Do m=1,k
           v(i,j) = v(i,j) + l(i,m)*u(m,j)
          End Do
        End Do
      End Do

      write (6,*)
      write (6,*) " verified a"
      write (6,*)

      Do i=1,n
        write (6,100) (v(i,j),j=1,n)
      End Do 

c------------
c Determinant
c------------

      Det = 1.0D0

      Do i=1,n
         Det = Det*d(i)
      End Do

      write (6,*) 
      write (6,101) Det

c-----
c Done
c-----

  100 Format (50(3x,f10.6))
  101 Format (" Determinant of A = ",f15.6)

      Stop
      End
