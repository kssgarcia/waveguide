      subroutine lobatto (m,Z)

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
c      Oxford University Press, 1998
c------------------------------------------------

c--------------------------------------------------
c Base points for (m+1)-point Lobatto interpolation
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Z(20)

c----------------------
c first and last points
c----------------------

      Z(1)   = -1.0D0
      Z(m+1) =  1.0D0
      
c--------------------
c intermediate points
c--------------------

      N = m-1

c--------------------
      If(N.eq.0) then
c--------------------

      Go to 99
 
c--------------------
      Else If(N.eq.1) then
c--------------------

      Z(2) = 0.0D0

c-------------------------
      Else If(N.eq.2) then
c-------------------------

      five = 5.0D0

      Z(2) = -1.0D0/Dsqrt(five)
      Z(3) = -Z(2)

c-------------------------
      Else If(N.eq.3) then
c-------------------------

      Z(2) = -0.65465 367
      Z(3) =  0.0D0
      Z(4) = -Z(2)

c-------------------------
      Else If(N.eq.4) then
c-------------------------

      Z(2) = -0.76505 532
      Z(3) = -0.28523 152
      Z(4) = -Z(3)
      Z(5) = -Z(2)

c-------------------------
      Else If(N.eq.5) then
c-------------------------

      Z(2) = -0.83022 390
      Z(3) = -0.46884 879 
      Z(4) =  0.0D0
      Z(5) = -Z(3)
      Z(6) = -Z(2)

c-------------------------
      Else If(N.eq.6) then
c-------------------------

      Z(2) = - 0.87174 015
      Z(3) = - 0.59170 018
      Z(4) = - 0.20929 922
      Z(5) = -Z(4)
      Z(6) = -Z(3)
      Z(7) = -Z(2)

c---------
      Else
c---------

      write (6,*) " lobatto: sorry, not yet implemented"
      stop

c-----------
      End If
c-----------
 
  99  Continue

      Return
      End
