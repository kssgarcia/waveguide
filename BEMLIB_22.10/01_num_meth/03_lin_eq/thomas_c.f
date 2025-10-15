      subroutine thomas_c 
     +
     +  (N   ! matrix size
     +  ,a   ! diagonal
     +  ,b   ! super-diagonal
     +  ,c   ! sub-diagonal
     +  ,s   ! right-hand side
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c          C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c------------------------------------------------

c------------------------------------------
c Compact Thomas algorithm for tridiagonal systems
c based on Gauss elimination
c
c Coefficient matrix:
c
c  | a1 b1  0   0  ...  0   0    0    |
c  | c2 a2  b2  0  ...  0   0    0    |
c  | 0  c3  a3  b3 ...  0   0    0    |
c  | ..............................   |
c  | 0  0   0   0  ... cn-1 an-1 bn-1 |
c  | 0  0   0   0  ...  0   cn   an   |
c
c------------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension a(1026),b(1026),c(1026),s(1026)
      Parameter (tol=0.000000001)

c----------
c reduction
c----------

      Do i=1,n-1
       i1 = i+1
       c(i1) = c(i1)/a(i)
       a(i1) = a(i1) - c(i1)*b(i)
       s(i1) = s(i1) - c(i1)*s(i)
      End Do

c------------------
c back substitution
c------------------

      s(N) = s(N)/a(N)

      Do i=Na,1,-1
        s(i)= (s(i)-b(i)*s(i+1))/a(i)
      End Do

c-----
c Done
c-----

 100  Format (1x,f15.10)

      Return
      End
