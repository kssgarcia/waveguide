      subroutine chol_c (l,a,n)

c=========================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c----------------------------------------
c This program accompanies the book:
c        C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press
c-----------------------------------------------

c---------------------------------
c Cholesky decomposition by column
c---------------------------------

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

c-------
c launch
c-------

      l(1,1) = Dsqrt(a(1,1))

      Do i=2,n
        l(i,1) = a(i,1)/l(1,1)
      End Do

      Do k=2,n

        sum = a(k,k)
        Do m=1,k-1
          sum= sum-l(k,m)*l(k,m)
        End Do
        l(k,k)=Dsqrt(sum)

        Do i=k+1,n
          sum = a(i,k)
          Do m=1,k-1
            sum = sum-l(i,m)*l(k,m)
          End Do
            l(i,k) = sum/l(k,k)
        End Do

      End Do

c-----
c done
c-----

      return
      end
