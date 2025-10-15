      subroutine gauss_log (NQ,Z,W)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c            C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c       Oxford University Press, 1998
c------------------------------------------------

c----------------------------------------
c Abscissae and weights for the Gauss quadrature
c for an integrand with a logarithmic singularity
c----------------------------------------
 
      Implicit Double Precision (a-h,o-z)

      Dimension Z(20),W(20)
 
c---
c set number of points
c---

      If(     NQ.ne.1
     +   .and.NQ.ne.2
     +   .and.NQ.ne.3
     +   .and.NQ.ne.4
     +   .and.NQ.ne.5
     +  ) then

        write (6,*) ' gauss_log: selected number of points '
        write (6,*) '            is not available; will take NQ=5'
        NQ=5 

      End If

c--------------------
      If(NQ.eq.1) then
c--------------------

      Z(1) = 0.25D0
      W(1) = 1.00D0
 
c--------------------
      Else If(NQ.eq.2) then
c--------------------

      Z(1) = 0.11200 88062 D0
      Z(2) = 0.60227 69081 D0

      W(1) = 0.71853 93190 D0
      W(2) = 0.28146 06809 D0

c-------------------------
      Else If(NQ.eq.3) then
c-------------------------

      Z(1) = 0.06389 07930 8 D0
      Z(2) = 0.36899 70637   D0
      Z(3) = 0.76688 03039   D0

      W(1) = 0.51340 45522   D0
      W(2) = 0.39198 00412   D0
      W(3) = 0.09461 54065 6 D0

c--------------------
      Else If(NQ.eq.4) then
c--------------------

      Z(1) = 0.04144 84801 9938322 D0
      Z(2) = 0.24527 49143 2060225 D0
      Z(3) = 0.55616 54535 6027584 D0
      Z(4) = 0.84898 23945 3298517 D0

      W(1) = 0.38346 40681 4513512 D0
      W(2) = 0.38687 53177 7476263 D0
      W(3) = 0.19043 51269 5014242 D0
      W(4) = 0.03922 54871 2995983 D0

c--------------------
      Else If(NQ.eq.5) then
c--------------------

      Z(1) = 0.02913 44721 51972 D0
      Z(2) = 0.17397 72133 208976287 D0
      Z(3) = 0.41170 25202 849020 D0
      Z(4) = 0.67731 41745 82820380 D0
      Z(5) = 0.89477 13610 3100828363 D0

      W(1) = 0.29789 34717 828944 D0
      W(2) = 0.34977 62265 1322418037 D0
      W(3) = 0.23448 82900 4405241888 D0
      W(4) = 0.09893 04595 16633146976 D0
      W(5) = 0.01891 15521 43195796489 D0

c-----------
      End If
c-----------

      Return
      End
