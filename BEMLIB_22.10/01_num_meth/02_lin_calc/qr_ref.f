	subroutine qr_reflex (a,n,q,r)

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
c             C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press, 1998
c------------------------------------------------

c--------------------------------------------------
c  Factorization of an arbitrary square matrix 'a'
c  into the product of an orthogonal matrix, q,
c  and an upper triangular matrix, r,
c  using successive Householder transformations.
c
c  Algorithm 2.7.1
c
c  SYMBOLS:
c  -------
c
c  a .... arbitrary square matrix
c  n .... size (rows/columns) of matrix a
c  q .... orthogonal matrix of size n
c  r .... upper triangular matrix of size n
c
c  del .. kronecker delta
c  w .... vector used in constructing q
c  u .... column vector from matrix r
c  c .... temporary vector matrix product
c  t .... storage matrix
c
c  eps... tolerance for skipping a reduction pass
c
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),q(128,128),r(128,128),del(128,128)
      Dimension w(129),u(128),c(128),t(128,128)

      Parameter (eps=0.0000001)

c------------------
c identity matrix
c------------------

      Do i=1,n
       Do j=1,n
         del(i,j) = 0.0D0
       End Do
       del(i,i) = 1.0D0
      End Do

c------------------------------------------------
c Copy the matrix a onto r, 
c which will become upper triangular
c
c Initialize the matrix q as the identity matrix.
c------------------------------------------------

      Do i=1,n
       Do j=1,n
         r(i,j) = a(i,j)
         q(i,j) = 0.0D0
       End Do
       q(i,i) = 1.0D0
      End Do

c-------------------------------------
c loop over the columns of the matrix r
c performing the Householder reduction
c-------------------------------------

      Do m=1,n-1          !  outer loop over columns of r

c---
c establish vector u as the ith column or r
c---

       Do j=m,n
         u(j)=r(j,m)
       End Do

c---
c compute the components of vector w
c---
 
        sum=0.0D0
 
        Do i=m,n
          sum=sum+u(i)**2
        End Do

        s=Dsqrt(sum)

        If(u(m).ge.0) then
          w(m) = Dsqrt(0.5D0*(1.0D0+u(m)/s))
        Else 
          w(m) = Dsqrt(0.5D0*(1.0D0-u(m)/s))
          s = -s
        End If

        If(abs(w(m)).gt.eps) then   !  skip if elements are small

        den = 2.0D0*s*w(m)

        Do i=m+1,n
          w(i)=u(i)/den
        End Do

        Do i=1,m-1
          w(i)=0.0D0
        End Do

c------------------------
c update matrices r and q
c------------------------

        Do j=m,n
          c(j)=0.0D0
          Do k=m,n
           c(j)=c(j)+w(k)*r(k,j)
          End Do
        End Do

        Do j=m,n
          Do k=m,n
           r(j,k)=r(j,k)-2.0D0*w(j)*c(k)
          End Do
        End Do

        Do j=1,n
          Do k=1,n
           t(j,k)=0.0D0
           Do l=1,n
            t(j,k)=t(j,k)+q(j,l)*(del(l,k)-2.0D0*w(l)*w(k))
           End Do
          End Do
        End Do

        Do j=1,n
         Do k=1,n
          q(j,k)=t(j,k)
         End Do
        End Do

       End If            !  end of small w(m) if block

c----------
      End Do              !  end of loop over columns of r
c----------

c-----
c Done
c-----

      Return
      End
