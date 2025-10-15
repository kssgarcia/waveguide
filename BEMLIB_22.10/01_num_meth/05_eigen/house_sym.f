      subroutine house_sym (n,a,b)

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
c            C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c----------------------------------------
c  Transforms a real symmetric matrix
c  to a real symmetric tridiagonal matrix
c  by Householder reflections.
c
c  Algorithm (5.7.19)
c
c  SYMBOLS:
c  --------
c
c  a .... real symmetric matrix
c  n .... size (rows/columns) of matrix a
c  b .... tri-diagonal matrix of size n
c
c  s .... sum of row elements excluding diagonal
c  r .... parameter for vector computation
c  w .... vector used in computation of updated matrix b
c  v .... vector used in computation of updated matrix b
c  z .... vector used in computation of updated matrix b
c
c  eps... tolerance
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(10,10),b(10,10),u(10),v(10),z(10)

      Parameter (eps=0.0000000001)

c----------------------
c copy matrix a 
c into working matrix b
c----------------------

      Do i=1,n
        Do j=1,n
         b(i,j)=a(i,j)
        End Do
      End Do

c-----------------------------------
c begin outer loop over sub-matrices
c from size n to 3
c-----------------------------------

      Do 95 i=1,n-2             !  outer loop over sub-matrices

c---
c compute sum of squares of elements
c on the working row,
c excluding the diagonal
c---

c---
c counterpart of equation (5.7.18)
c with i in place of 1
c---

      sum = 0.0D0

      Do j=i+1,n
        sum = sum+b(i,j)**2
      End Do

      test = (n-i)*eps
      If(sum.lt.test) Go to 95   !  element of this row are
                                 !  virtually zero; skip this row
  
c---
c proceed with the reflection
c---

      check = b(i,i+1)

      prosimo = check/Dabs(check)
      s = prosimo * Dsqrt(sum)

      den  = s*(s+b(i,i+1))
      den2 = 2.0D0*den

c---
c compute the vector u
c---
  
      Do j=i,n

        If(j.eq.i) then
          u(j) = 0.0D0
        Else If(j.eq.i+1) then
          u(j) = b(i,i+1)+s
        Else
         u(j) = b(i,j)
        End If

      End Do

c---
c compute the vector z
c---

      Do j=i,n
        sum = 0.0D0
        Do k=i+1,n
         sum = sum+b(j,k)*u(k)
        End Do 
        z(j) = sum/den
      End Do

c---
c compute the vector v
c---

      Do j=i,n
         sum = 0.0D0
         Do k=i+1,n
          sum = sum+u(k)*z(k)
         End Do
         v(j) = z(j)-sum*u(j)/den2
      End Do  

c---
c compute the components of the upper part of
c the lower right block of the updated matrix b
c---

       Do j=i,n
         Do k=i,n
           b(j,k) = b(j,k)-u(j)*v(k)-u(k)*v(j)
         End Do
       End Do

c---
c set the components of the lower part of
c the lower right block of the updated matrix b
c---

       Do j=i,n
         Do k=i,n
           b(k,j) =  b(j,k) 
         End Do
       End Do

  95  Continue           !  end loop over sub matrices
  
c-----
c Done
c-----

      Return
      End
