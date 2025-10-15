      subroutine cramer_33 
     +
     +  (A11,A12,A13
     +  ,A21,A22,A23
     +  ,A31,A32,A33
     +  ,B1,B2,B3
     +
     +  ,X1,X2,X3
     +  )

c===========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c===========================================

c------------------------------------------------
c This program accompanies the book:
c
c            C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c          Oxford University Press
c------------------------------------------------

c--------------------------------------------
c Solves a 3x3 real system by Cramer's rule
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

c---
c Cramer's rule
c---

      Det =  A11*( A22*A33-A23*A32 )
     +     - A12*( A21*A33-A23*A31 )
     +     + A13*( A21*A32-A22*A31 )

      Det1 =  B1*( A22*A33-A23*A32 )
     +     - A12*(  B2*A33-A23*B3  )
     +     + A13*(  B2*A32-A22*B3  )

      Det2 = A11*( B2 *A33-A23*B3  )
     +     -  B1*( A21*A33-A23*A31 )
     +     + A13*( A21* B3-B2 *A31 )

      Det3 = A11*( A22* B3-A32* B2 )
     +     - A12*( A21* B3-A31* B2 )
     +     +  B1*( A21*A32-A22*A31 )

      X1 = Det1/Det
      X2 = Det2/Det
      X3 = Det3/Det

c----------
c residuals
c---
c
c      res1 = B1 - A11*X1-A12*X2-A13*X3
c      res2 = B2 - A21*X1-A22*X2-A23*X3
c      res3 = B3 - A31*X1-A32*X2-A33*X3
c      write (6,*), res1,res2,res3
c----------

      return
      end
