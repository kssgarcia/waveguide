      subroutine chol_r (l,a,n)

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c==========================================

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c------------------------------------------------

c------------------------------
c Cholesky decomposition by row
c Algorithm 2.6.4
c------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision l(10,10)
      Dimension a(10,10)

c-----------
c initialize
c----------

      Do i=1,n
       Do j=1,n
        l(i,j) = 0.0D0
       End Do
      End Do

c---------
c Cholesky
c---------

      l(1,1) = Dsqrt(a(1,1))

      Do i=2,n

         Do j=1,i-1
           sum = a(i,j)
           Do m=1,j-1
              sum=sum-l(i,m)*l(j,m)
           End Do
           l(i,j)=sum/l(j,j)
         End Do

         sum = a(i,i)
         Do m=1,i-1
           sum = sum-l(i,m)*l(i,m)
         End Do
         l(i,i) = Dsqrt(sum)

      End Do

c-----
c done
c-----

      return
      end
