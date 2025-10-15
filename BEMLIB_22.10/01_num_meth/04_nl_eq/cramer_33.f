      subroutine cramer33 
     +
     + (A11,A12,A13
     + ,A21,A22,A23
     + ,A31,A32,A33
     + ,B1,B2,B3
     + ,X1,X2,X3
     + )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c Solution of a 3x3 real system by Cramer's rule
c
c This may appear to be a silly program,
c but keep it; some time, it may save
c you a lot of time
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

c--------------
c Cramer's rule
c--------------

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

c---
c verify
c---
c      Check1 = B1 - A11*X1-A12*X2-A13*X3
c      Check2 = B2 - A21*X1-A22*X2-A23*X3
c      Check3 = B3 - A31*X1-A32*X2-A33*X3
c      Print *, Check1,Check2,Check3
c---

      Return
      End
