      subroutine cramer_33c
     +
     +   (AR,AI
     +   ,BR,BI
     +   ,XR,XI
     +   )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c             C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press
c------------------------------------------------

c---------------------------------------
c  Solution of a complex 3x3 system
c  by Cramer's rule
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

      Dimension AR(3,3),AI(3,3)

      Dimension BR(3),BI(3)
      Dimension XR(3),XI(3)

      Dimension detR(3),detI(3)
      Dimension savR(3),savI(3)

      Parameter (tol=0.0000001D0)

c--------
c display
c--------

      Go to 333

      write (6,*)
      write (6,101) AR(1,1),AR(1,2),AR(1,3),BR(1)
      write (6,101) AR(2,1),AR(2,2),AR(2,3),BR(2)
      write (6,101) AR(3,1),AR(3,2),AR(3,3),BR(3)
      write (6,*)
      write (6,101) AI(1,1),AI(1,2),AI(1,3),BI(1)
      write (6,101) AI(2,1),AI(2,2),AI(2,3),BI(2)
      write (6,101) AI(3,1),AI(3,2),AI(3,3),BI(3)

 333  Continue

c------------------------
c Compute the determinant
c------------------------

      call det_33c
     +
     +   (AR,AI
     +   ,detR0,detI0
     +   )

c-----------------
c loop over minors
c-----------------

      Do k=1,3

        Do j=1,3
         savR(j) = AR(j,k)
         savI(j) = AI(j,k)
         AR(j,k) = BR(j)
         AI(j,k) = BI(j)
        End Do

        call det_33c
     +
     +   (AR,AI
     +   ,detR(k),DetI(k)
     +   )

c       restore:

        Do j=1,3
         AR(j,k) = savR(j)
         AI(j,k) = savI(j)
        End Do

      End Do

c--------------
c Cramer's rule
c--------------

      den =  detR0*detR0 + detI0*detI0
      dtR =  detR0/den 
      dtI = -detI0/Den 

      Do k=1,3
       XR(k) = DetR(k)*DtR - DetI(k)*DtI
       XI(k) = DetR(k)*DtI + DetI(k)*DtR
      End Do

c-------
c verify
c-------

      Do k=1,3
       VR = 0.0D0
       VI = 0.0D0
       Do j=1,3
        VR = VR + AR(k,j)*XR(J)-AI(k,j)*XI(j)
        VI = VI + AR(k,j)*XI(J)+AI(k,j)*XR(j)
       End Do
       VR = VR - BR(k)
       VI = VI - BI(k)

       If(Dabs(VR).gt.tol.or.Dabs(VI).gt.tol) then
        write (6,*) " cramer_33c: failed",Dabs(VR),Dabs(VI)
       End If

c     write (6,*) " cramer_33c: succesful"

      End Do
      
c-----
c Done
c-----

  100 Format (1x,I3,10(1x,f10.5))
  101 Format (10(1x,f10.5))

      Return
      End
