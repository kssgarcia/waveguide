      subroutine qr_rot (a,n,q,r)

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
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c       Oxford University Press, 1998
c------------------------------------------------

c------------------------------------------------
c  Factorization of the arbitrary square matrix A
c  into the product of an orthogonal matrix, Q,
c  and the upper triangular matrix, R,
c  by successive Givens transformations.
c
c  Algorithm 2.7.2
c
c  SYMBOLS:
c  -------
c
c  a .... arbitrary square matrix
c  n .... size (rows/columns) of matrix a
c  q .... orthogonal matrix of size n
c  r .... upper triangular matrix of size n
c
c  s .... sin component of orthogonal transformation
c  c .... cos component of orthogonal transformation
c  p .... parameter computed from elements of r
c
c  eps... tolerance for skipping a reduction
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),q(128,128),r(128,128)
      Parameter (eps=0.0000001)

c------------------------------------------
c Copy matrix a onto the working matrix r,
c which will become upper triangular
c
c Initialize the orthogonal matrix as
c the identity matrix
c------------------------------------------

      Do i=1,n
       Do j=1,n
        r(i,j)=a(i,j)
        q(i,j)=0.0D0
       End Do
       q(i,i)=1.0D0
      End Do

c-----------------------------------------------------
c  loop over columns and rows of matrix r successively
c  zeroing elements via Givens rotations, and perform
c  transformations on the matrix q.
c-----------------------------------------------------

      Do j=1,n-1                 !  working column of matrix r

       Do i=j+1,n                !  working row of matrix r
    
         If(Dabs(r(i,j)).gt.eps) then  ! skip if element is small  

          p = Dsqrt(r(j,j)**2+r(i,j)**2)
          c = r(j,j)/p
          s = r(i,j)/p

          Do k=j,n
           t1 =  c*r(j,k)+s*r(i,k)
           t2 = -s*r(j,k)+c*r(i,k)
           r(j,k) = t1
           r(i,k) = t2
          End Do

          Do k=1,n
           t1 =  c*q(j,k)+s*q(i,k)
           t2 = -s*q(j,k)+c*q(i,k)
           q(j,k) = t1
           q(i,k) = t2
          End Do

          End If

       End Do                   !  end of row loop
      End Do                    !  end of column loop

c---------------------------------------
c New matrix q is the transpose of old q
c---------------------------------------

      Do i=1,n-1
        Do j=i+1,n
          temp = q(i,j)
          q(i,j) = q(j,i)
          q(j,i) = temp
        End Do
      End Do

c-----
c Done
c-----

      Return
      End
