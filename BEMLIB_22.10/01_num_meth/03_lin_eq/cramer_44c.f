      subroutine cramer_44c
     +
     +   (AR,AI
     +   ,BR,BI
     +   ,XR,XI
     +   )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c----------------------------------------------------
c This program accompanies the book:
c
c             C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c       Oxford University Press
c----------------------------------------------------

c---------------------------------------
c  Solution of a 4x4 complex system
c  using Cramer's rule
c
c  LEGEND:
c  ------
c
c  AR:  real part of the matrix
c  AI:  imaginary part of the matrix
c
c  BR:  real part of the rhs
c  BI:  imaginary part of the rhs
c
c  XR:  real part of the solution
c  XI:  imaginary part of the solution
c
c  detR:  real part of the determinant
c  detI:  imaginary part of the determinant
c------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension AR(4,4),AI(4,4)

      Dimension BR(4),BI(4)
      Dimension XR(4),XI(4)

      Dimension detR(4),detI(4)
      Dimension savR(4),savI(4)

      Parameter (tol=0.000000001)

c--------
c display
c--------

      Go to 333

      write (6,*)
      write (6,101) AR(1,1),AR(1,2),AR(1,3),AR(1,4),BR(1)
      write (6,101) AR(2,1),AR(2,2),AR(2,3),AR(2,4),BR(2)
      write (6,101) AR(3,1),AR(3,2),AR(3,3),AR(3,4),BR(3)
      write (6,101) AR(4,1),AR(4,2),AR(4,3),AR(4,4),BR(4)
      write (6,*)
      write (6,101) AI(1,1),AI(1,2),AI(1,3),AI(1,4),BI(1)
      write (6,101) AI(2,1),AI(2,2),AI(2,3),AI(2,4),BI(2)
      write (6,101) AI(3,1),AI(3,2),AI(3,3),AI(3,4),BI(3)
      write (6,101) AI(4,1),AI(4,2),AI(4,3),AI(4,4),BI(4)

 333  Continue

c-----------
c Compute the determinant
c-----------

      call det_44c
     +
     +   (AR,AI
     +   ,detR0,detI0
     +   )

c---
c loop over minors
c---

      Do k=1,4

        Do j=1,4
         savR(j) = AR(j,k)
         savI(j) = AI(j,k)
         AR(j,k) = BR(j)
         AI(j,k) = BI(j)
        End Do

        call det_44c
     +
     +   (AR,AI
     +   ,detR(k),detI(k)
     +   )

        Do j=1,4
         AR(j,k) = savR(j)
         AI(j,k) = savI(j)
        End Do

      End Do

c-----
c Cramer's rule
c-----

      den =  detR0*detR0 + detI0*detI0
      DtR =  detR0/den 
      DtI = -detI0/den 

      Do k=1,4
       XR(k) = detR(k)*DtR - detI(k)*DtI
       XI(k) = detR(k)*DtI + detI(k)*DtR
      End Do

c-----
c verify
c-----

      Do k=1,4

       VR = 0.0D0
       VI = 0.0D0

       Do j=1,4
        VR = VR + AR(k,j)*XR(j)-AI(k,j)*XI(j)
        VI = VI + AR(k,j)*XI(j)+AI(k,j)*XR(j)
       End Do
       VR = VR - BR(k)
       VI = VI - BI(k)

       If(Dabs(VR).gt.tol.or.Dabs(VI).gt.tol) then
        write (6,*) " cramer_44c: failed"
       End If

      End Do
      
c-----
c Done
c-----

  100 Format (1x,I3,10(1x,f10.5))
  101 Format (10(1x,f10.5))

      Return
      End
