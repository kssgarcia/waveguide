      subroutine Det_33c 
     + 
     +     (AR,AI
     +     ,DetR,DetI
     +     )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c "Numerical Computation in Science and Engineering"
c        Oxford University Press
c------------------------------------------------

c---------------------------------------
c  DETERMINANT OF A 3x3 COMPLEX MATRIX 
c
c  SYMBOLS:
c  -------
c
c  AR:	real part of the matrix
c  AI:	imaginary part of the matrix
c
c  DetR:  real part of the determinant
c  DetI:  imaginary part of the determinant
c------------------------------------------

      Implicit Double Precision (A-H,O-Z)
      Dimension AR(3,3),AI(3,3)

c-----------------
c A lot of algebra
c-----------------

      SMR1 = AR(2,2)*AR(3,3)-AI(2,2)*AI(3,3)
      SMI1 = AR(2,2)*AI(3,3)+AI(2,2)*AR(3,3)
      SMR2 = AR(2,3)*AR(3,2)-AI(2,3)*AI(3,2)
      SMI2 = AR(2,3)*AI(3,2)+AI(2,3)*AR(3,2)
      SMR  = SMR1-SMR2
      SMI  = SMI1-SMI2
      DTR1 = AR(1,1)*SMR - AI(1,1)*SMI
      DTI1 = AR(1,1)*SMI + AI(1,1)*SMR

      SMR1 = AR(2,1)*AR(3,3)-AI(2,1)*AI(3,3)
      SMI1 = AR(2,1)*AI(3,3)+AI(2,1)*AR(3,3)
      SMR2 = AR(2,3)*AR(3,1)-AI(2,3)*AI(3,1)
      SMI2 = AR(2,3)*AI(3,1)+AI(2,3)*AR(3,1)
      SMR  = SMR1-SMR2
      SMI  = SMI1-SMI2
      DTR2 = AR(1,2)*SMR - AI(1,2)*SMI
      DTI2 = AR(1,2)*SMI + AI(1,2)*SMR

      SMR1 = AR(2,1)*AR(3,2)-AI(2,1)*AI(3,2)
      SMI1 = AR(2,1)*AI(3,2)+AI(2,1)*AR(3,2)
      SMR2 = AR(2,2)*AR(3,1)-AI(2,2)*AI(3,1)
      SMI2 = AR(2,2)*AI(3,1)+AI(2,2)*AR(3,1)
      SMR  = SMR1-SMR2
      SMI  = SMI1-SMI2
      DTR3 = AR(1,3)*SMR - AI(1,3)*SMI
      DTI3 = AR(1,3)*SMI + AI(1,3)*SMR

      DetR = DTR1-DTR2+DTR3
      DetI = DTI1-DTI2+DTI3

c-----
c Done
c-----

      Return
      End
