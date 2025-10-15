      subroutine Det_44c 
     +
     +   (AR,AI
     +   ,DetR,DetI
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
c C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c Oxford University Press, 1998
c------------------------------------------------

c---------------------------------------
c
c  DETERMINANT OF A 4x4 COMPLEX MATRIX
c
c  AR:	real part of the matrix
c  AI:	imaginary part of the matrix
c
c  DetR:  real part of the determinant
c  DetI:  imaginary part of the determinant
c
c------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      Dimension AR(4,4),AI(4,4)

      Dimension BR(3,3),BI(3,3)

c-----------
c initialize
c-----------

      DetR = 0.0D0
      DetI = 0.0D0

c---
c Do it the long way by the 
c Laplace expansion
c---

      elmR = AR(1,1)
      elmI = AI(1,1)

      BR(1,1) = AR(2,2)
      BR(1,2) = AR(2,3)
      BR(1,3) = AR(2,4)
      BR(2,1) = AR(3,2)
      BR(2,2) = AR(3,3)
      BR(2,3) = AR(3,4)
      BR(3,1) = AR(4,2)
      BR(3,2) = AR(4,3)
      BR(3,3) = AR(4,4)

      BI(1,1) = AI(2,2)
      BI(1,2) = AI(2,3)
      BI(1,3) = AI(2,4)
      BI(2,1) = AI(3,2)
      BI(2,2) = AI(3,3)
      BI(2,3) = AI(3,4)
      BI(3,1) = AI(4,2)
      BI(3,2) = AI(4,3)
      BI(3,3) = AI(4,4)

      call Det_33c (BR,BI,DTR,DtI)

      DetR = DetR + elmR*DtR - elmI*DtI
      DetI = DetI + elmR*DtI + elmI*DtR

c---
      elmR  = - AR(2,1)
      elmI  = - AI(2,1)

      BR(1,1) = AR(1,2)
      BR(1,2) = AR(1,3)
      BR(1,3) = AR(1,4)
      BR(2,1) = AR(3,2)
      BR(2,2) = AR(3,3)
      BR(2,3) = AR(3,4)
      BR(3,1) = AR(4,2)
      BR(3,2) = AR(4,3)
      BR(3,3) = AR(4,4)

      BI(1,1) = AI(1,2)
      BI(1,2) = AI(1,3)
      BI(1,3) = AI(1,4)
      BI(2,1) = AI(3,2)
      BI(2,2) = AI(3,3)
      BI(2,3) = AI(3,4)
      BI(3,1) = AI(4,2)
      BI(3,2) = AI(4,3)
      BI(3,3) = AI(4,4)

      call Det_33c (BR,BI,DTR,DtI)

      DetR = DetR + elmR*DtR - elmI*DtI
      DetI = DetI + elmR*DtI + elmI*DtR

c---

      elmR = AR(3,1)
      elmI = AI(3,1)

      BR(1,1) = AR(1,2)
      BR(1,2) = AR(1,3)
      BR(1,3) = AR(1,4)
      BR(2,1) = AR(2,2)
      BR(2,2) = AR(2,3)
      BR(2,3) = AR(2,4)
      BR(3,1) = AR(4,2)
      BR(3,2) = AR(4,3)
      BR(3,3) = AR(4,4)

      BI(1,1) = AI(1,2)
      BI(1,2) = AI(1,3)
      BI(1,3) = AI(1,4)
      BI(2,1) = AI(2,2)
      BI(2,2) = AI(2,3)
      BI(2,3) = AI(2,4)
      BI(3,1) = AI(4,2)
      BI(3,2) = AI(4,3)
      BI(3,3) = AI(4,4)

      call Det_33c (BR,BI,DTR,DtI)

      DetR = DetR + elmR*DtR - elmI*DtI
      DetI = DetI + elmR*DtI + elmI*DtR

c---
      elmR = - AR(4,1)
      elmI = - AI(4,1)

      BR(1,1) = AR(1,2)
      BR(1,2) = AR(1,3)
      BR(1,3) = AR(1,4)
      BR(2,1) = AR(2,2)
      BR(2,2) = AR(2,3)
      BR(2,3) = AR(2,4)
      BR(3,1) = AR(3,2)
      BR(3,2) = AR(3,3)
      BR(3,3) = AR(3,4)

      BI(1,1) = AI(1,2)
      BI(1,2) = AI(1,3)
      BI(1,3) = AI(1,4)
      BI(2,1) = AI(2,2)
      BI(2,2) = AI(2,3)
      BI(2,3) = AI(2,4)
      BI(3,1) = AI(3,2)
      BI(3,2) = AI(3,3)
      BI(3,3) = AI(3,4)

      call Det_33c (BR,BI,DTR,DtI)

      DetR = DetR + elmR*DtR - elmI*DtI
      DetI = DetI + elmR*DtI + elmI*DtR

c-----
c Done
c-----

      Return
      End

c=================================================

      subroutine Det_33c (AR,AI,DetR,DetI)

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
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press, 1998
c
c------------------------------------------------

c------
c
c  DETERMINANT OF A 3x3 COMPLEX MATRIX A
c
c  AR:	real part of the matrix
c  AI:	imaginary part of the matrix
c
c  DetR:	real part of the determinant
c  DetI:	imaginary part of the determinant
c
c------

      Implicit Double Precision (A-H,O-Z)

      Dimension AR(3,3),AI(3,3)

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
