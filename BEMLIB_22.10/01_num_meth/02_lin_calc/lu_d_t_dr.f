      program lu_d_t_dr

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

c---------------------------------------
c Driver for the Doolittle decomposition
c of a tridiagonal matrix
c---------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension t(10,10),u(10,10),v(10,10)
      Dimension a(10),b(10),c(10)

      Double precision l(10,10)

c--------------------------
c read the diagonal vectors
c--------------------------

      open (unit=8,file='matrix_t.dat')
        read (8,*) n
        read (8,*) (a(i), i=1,n)
        read (8,*) (b(i), i=1,n-1)
        read (8,*) (c(i), i=2,n)
      close (8)

c---
c generate the tridiagonal matrix
c---

      Do i=1,n
       Do j=1,n
        t(i,j) = 0.0D0
       End Do
      End Do

      t(1,1) = a(1)
      t(1,2) = b(1)

      Do i=2,n
        t(i,i)   = a(i)
        t(i,i+1) = b(i)
        t(i,i-1) = c(i)
      End Do 

c---
c print
c---

      write (6,*)
      write (6,*) " Tridiagonal matrix"
      write (6,*)

      Do i=1,n
        write (6,1000) (t(i,j),j=1,n)
      End Do 
	
c----------
c Decompose
c---------

      call lu_d_t (n,c,a,b,l,u)

c---
c Print the three matrices
c---

      write (6,*)
      write (6,*) " Matrix L"
      write (6,*)

      Do i=1,n
        write (6,1000) (l(i,j),j=1,n)
      End Do 

      write (6,*)
      write (6,*) " Matrix U"
      write (6,*)

      Do i=1,n
        write(6,1000) (u(i,j),j=1,n)
      End Do 

c---
c verify the decomposition
c---

      Do i=1,n
        Do j=1,n
          v(i,j) = 0.0D0
          Do m=1,n
            v(i,j) = v(i,j)+l(i,m)*u(m,j)
          End Do
        End Do
      End Do 

      write (6,*)
      write (6,*) " Verified matrix T"
      write (6,*)

      Do i=1,n
        write(6,1000) (v(i,j),j=1,n)
      End Do 

c---
c Determinant
c---

      Det = 1.0D0
      Do i=1,n
        Det = Det*u(i,i)
      End Do
      write (6,*) " Determinant = ",Det

c-----
c Done
c-----

1000  format(10(3x,f10.7))

      Stop
      End
