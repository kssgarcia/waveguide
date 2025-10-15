      subroutine gauss_her (NQ,Z,W)

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c===================================================
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c===================================================

c========================================================
c Abscissae and weights for the Gauss-Hermitte quadrature
c applicable to a Gaussian-like function
c
c The weights add to sqrt(pi) for any NQ
c========================================================
 
      Implicit Double Precision (a-h,o-z) 

      Dimension Z(20),W(20)
 
c---
c set the number of points
c---

      If(     NQ.ne.1
     +   .and.NQ.ne.2
     +   .and.NQ.ne.3
     +   .and.NQ.ne.4
     +   .and.NQ.ne.5
     +  ) then
        write (6,*) ' gauss_her: Number of points is not available'
        write (6,*) '            Will take NQ=5'
        NQ = 5
      End If

c----------
c constants
c----------

      pi = 3.14159265358D0
      srpi = Dsqrt(pi)
 
c--------------------
      if(NQ.eq.1) then
c--------------------

      Z(1) = 0.0D0
      W(1) = srpi

c--------------------------
      else if(NQ.eq.2) then
c--------------------------

      Z(1) =-1.0D0/Dsqrt(2.0D0)
      Z(2) =-Z(1)

      W(1) = 0.5D0*srpi
      W(2) = W(1)

c--------------------------
      else if(NQ.eq.3) then
c--------------------------

      Z(1) =-Dsqrt(1.5D0)
      Z(2) = 0.0D0
      Z(3) =-Z(1)

      W(1) = 0.2954089752D0
      W(2) = 1.181635901D0

      fc = srpi*6.0D0*(2.0D0**(NQ+1))

      Do i=1,2
       W(i) = fc/(16.0D0*Z(i)**4-48.0D0*Z(i)**2+12.0D0)**2
      End Do

      W(3) = W(1)

c--------------------------
      else if(NQ.eq.4) then
c--------------------------

      Z(1) = -1.65068012388578D0
      Z(2) = -0.52464762327529D0
      Z(3) = -Z(2)
      Z(4) = -Z(1)

      W(1) = 0.08131283545 D0
      W(2) = 0.8049140900 D0

      fc = srpi*24.0D0*(2.0D0**(NQ+1))

      Do i=1,2
       W(i) = fc/(32.0D0*Z(i)**5-160.0D0*Z(i)**3+120.0D0*Z(i))**2
      End Do

      W(3) = W(2)
      W(4) = W(1)

c--------------------------
      else if(NQ.eq.5) then
c--------------------------

      Z(1) =-2.020182870D0
      Z(2) =-0.9585724646D0 
      Z(3) =  0.0D0
      Z(4) =-Z(2)
      Z(5) =-Z(1)

      W(1) = 0.01995324206D0
      W(2) = 0.3936193232D0
      W(3) = 0.9453087205D0

      fc = srpi*120.0D0*(2.0D0**(NQ+1))

      Do i=1,3
       W(i) = fc/(64.0D0*Z(i)**6-480.0D0*Z(i)**4+720.0D0*Z(i)**2
     +           -120.0D0)**2
      End Do

      W(4) = W(2)
      W(5) = W(1)

c-----------
      end if
c-----------

      Return
      End
