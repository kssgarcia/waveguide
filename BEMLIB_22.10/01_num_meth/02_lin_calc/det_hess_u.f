      subroutine det_hess_u (n,a,det)

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

c----------------------------------------
c Computation of the determinant of
c an upper Hessenberg matrix
c according to algorithm (2.4.11)
c
c SYMBOLS:
c --------
c
c  a .... upper hessenberg matrix
c  n .... size (rows/columns) of matrix a
c  det .. determinant of matrix a
c	
c  u .... upper triangular matrix constructed from a
c  d .... diagonal matrix constructed from a
c----------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension a(10,10),u(10,10),d(0:10),p(0:10)

c---
c generate the diagonal matrix and u
c---

      Do i=1,n-1
        d(i) = - a(i+1,i)
      End Do

      d(n) = 1.0D0

      Do i=1,n
        Do j=i,n
         u(i,j) = a(i,j)/d(j)
        End Do
      End Do

c---
c  compute the determinant 
c---

      p(0) = 1.0D0
      p(1) = u(1,1)

      Do i=2,n
        p(i) = 0.0D0
        Do j=1,i
          p(i) = p(i)+u(j,i)*p(j-1)
        End Do
      End Do

c---
c  the determinant of the matrix a is the product of d(n)
c  and p(n)
c---

      det = p(n)

      Do i=1,n-1
        det = det*d(i)
      End Do

c-----
c Done
c-----

      Return
      End
