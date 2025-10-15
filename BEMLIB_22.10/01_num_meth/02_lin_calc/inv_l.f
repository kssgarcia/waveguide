      subroutine inv_l (n,a,b)

c-----------------------------------------
c Copyright 1999 by C. Pozrikidis
c       All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c Numerical Computation in Science and Engineering
c        Oxford University Press
c------------------------------------------------

c---------------------------------------------
c Inverse of an nxn lower triangular matrix: a
c Inverse is matrix: b
c---------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension a(128,128),b(128,128)

c-------
c launch
c-------

      Do i=1,N
       b(i,i) = 1.0D0/a(i,i)
      End Do

      Do j=1,N-1
        Do i=j+1,N
          sum = 0.0D0
          Do m=j,i-1
            sum = sum +a(i,m)*b(m,j)
          End Do
          b(i,j) = - sum/a(i,i)
        End Do
      End Do

c-----
c Done
c-----

      return
      end
