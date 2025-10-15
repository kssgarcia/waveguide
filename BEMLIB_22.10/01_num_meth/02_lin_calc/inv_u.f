      subroutine inv_u (n,a,b)

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
c         C. Pozrikidis
c Numerical Computation in Science and Engineering
c      Oxford University Press, 1998
c------------------------------------------------

c------------------------------------------------
c Computation of the inverse of nxn 
c upper triangular matrix: a
c
c Inverse is matrix: b
c
c Algorithm (2.5.8)
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(100,100),b(100,100)

c-------
c launch
c-------

      Do i=1,N
       b(i,i) = 1.0D0/a(i,i)
      End Do

      Do j=2,N
        Do i=j-1,1,-1
          sum = 0.0D0
          Do m=i+1,j
            sum = sum + a(i,m)*b(m,j)
          End Do
          B(i,j) = - sum/a(i,i)
        End Do
      End Do

c-----
c Done
c-----

      Return
      End
