      subroutine crout (n,a,l,u)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c          C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c     Oxford University Press, 1998
c------------------------------------------------

c---------------------
c  Crout decomposition
c
c  Algorithm 2.6.2
c---------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision l(100,100)
      Dimension        u(100,100),a(100,100)

c---------------------------
c first row and first column
c---------------------------

      Do i=1,n
        u(i,i) = 1.0D0
        l(i,1) = a(i,1)
        u(1,i) = a(1,i)/l(1,1)
      End Do

c-------
c others
c-------

      Do k=2,n

        Do i=k,n
         sum=0.0D0
         Do m=1,k-1
           sum = sum + l(i,m)*u(m,k)
         End Do
         l(i,k) = a(i,k) - sum
        End Do

        Do j=k+1,n
          sum = 0.0D0
          Do m=1,k-1
            sum = sum + l(k,m)*u(m,j)
          End Do
          u(k,j) = (a(k,j)-sum)/l(k,k)
        End Do

      End Do

c-----
c Done
c-----

      Return
      End
