      subroutine cramer_22
     +
     +   (A11,A12
     +   ,A21,A22
     +   ,B1,B2
     +   ,X1,X2
     +   )

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
c                C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press, 1998
c------------------------------------------------

c-----------------------------------------
c Solve a real 2x2 system by Cramer's rule
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

c---------------
c Cramer's rule:
c---------------

      Det =  A11*A22-A21*A12

      Det1 =   B1*A22 - B2*A12
      Det2 = - B1*A21 + B2*A11

      X1 = Det1/Det
      X2 = Det2/Det

c----------
c residuals
c---
c
c      res1 = B1 - A11*X1-A12*X2
c      res2 = B2 - A21*X1-A22*X2
c      write (6,*), res1,res2
c----------

c-----
c Done
c-----

      Return
      End
