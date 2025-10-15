      subroutine radau (NQ,Z,W)

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press
c------------------------------------------------

c------------------------------------------------
c Base points and weights for the Radau
c quadrature with NQ points
c
c Integral = sum_{i=1}^{NQ} f_i w_i
c
c The first base point is set at -1
c The intermediate NQ-1 points are located
c at the zeros of the (NQ-1)-degree Radau polynomial
c
c This table contains values for
c
c  NQ = 1,2,3,4,5,6
c
c  Default value is 6
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Z(20),W(20)

c------
c check
c------

      if(     NQ.ne.1
     +   .and.NQ.ne.2
     +   .and.NQ.ne.3
     +   .and.NQ.ne.4
     +   .and.NQ.ne.5
     +   .and.NQ.ne.6
     +  ) then

        write (6,*) "radau: selected number of points"
        write (6,*) "       is not available; will take NQ=6"
        NQ=5

      end if

      Z(1) =-1.0D0;
      W(1) = 2.0D0/NQ**2;

      fc = 0.5D0*W(1)

c--------------------
      if(NQ.eq.2) then
c--------------------

      Z(2) = 1.0D0/3.0D0
      W(2) = 3.0D0/2.0D0

c-------------------------
      else if(NQ.eq.3) then
c-------------------------

      Z(2) = -0.28989794855664 D0
      Z(3) =  0.68989794855664 D0

      Do i=2,3
       W(i) = fc*(1-Z(i))/( 0.5*(3.0D0*Z(i)**2-1.0D0) )**2
      End Do

c-------------------------
      else if(NQ.eq.4) then
c-------------------------

      Z(2) = -0.57531892352169 D0
      Z(3) =  0.18106627111853 D0
      Z(4) =  0.82282408097459 D0

      Do i=2,4
       W(i) = fc*(1-Z(i))/( 0.5*(5.0D0*Z(i)**3-3.0D0*Z(i)) )**2
      End Do

c-------------------------
      else if(NQ.eq.5) then
c-------------------------

      Z(2) = -0.72048027131244 D0
      Z(3) = -0.16718086473783 D0
      Z(4) =  0.44631397272375 D0
      Z(5) =  0.88579160777096 D0

      Do i=2,5
      W(i) = fc*(1-Z(i))
     +    /( (35.0D0*Z(i)**4-30.0D0*Z(i)**2+3.0D0)/8.0 )**2
      End Do

c-------------------------
      else if(NQ.eq.6) then
c-------------------------

      Z(2) = -0.80292982840235 D0
      Z(3) = -0.39092854670727 D0
      Z(4) =  0.12405037950523 D0
      Z(5) =  0.60397316425278 D0
      Z(6) =  0.92038028589706 D0

      Do i=2,6
      W(i) = fc*(1-Z(i))
     + /( (63.0D0*Z(i)**5-70.0D0*Z(i)**3+15.0D0*Z(i))/8.0 )**2
      End Do

c-----------
      end if
c-----------
 
      return
      end
