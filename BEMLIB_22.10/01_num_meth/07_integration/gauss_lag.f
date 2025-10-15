      subroutine gauss_lag (NQ,Z,W)

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c This program accompanies the book:
c           C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c        Oxford University Press
c------------------------------------------------

c------------------------------
c Abscissae and weights for the
c Gauss--Laguerre quadrature
c-------------------------------
 
      Implicit Double Precision (a-h,o-z) 

      Dimension Z(20),W(20)
 
c---
c set the number of points
c---

      if(     NQ.ne.1
     +   .and.NQ.ne.2
     +   .and.NQ.ne.3
     +   .and.NQ.ne.4
     +   .and.NQ.ne.5
     +  ) then
        write (6,*) ' gauss_lag: Number of points is not available'
        write (6,*) '            Will take NQ=5'
        NQ = 5
      end if

c---
c prepare
c---

      fc = 1.0D0

      Do i=1,NQ
       fc = fc*i*i;
      End Do
 
c--------------------
      if(NQ.eq.1) then
c--------------------

      Z(1) = 0.0D0
      W(1) = 1.0D0

c--------------------
      else if(NQ.eq.2) then
c--------------------

      Z(1) = 0.58578643762690 D0
      Z(2) = 3.41421356237309 D0

      W(1) = 0.8535533905 D0
      W(2) = 0.1464466094 D0

      Do i=1,2
       W(i) = fc*Z(i)/( -Z(i)**3+9.0D0*Z(i)**2-18.0D0*Z(i)+6.0D0 )**2
      End Do

c-------------------------
      else if(NQ.eq.3) then
c-------------------------

      Z(1) = 0.41577455678348 D0
      Z(2) = 2.29428036027904 D0
      Z(3) = 6.28994508293748 D0

      W(1) = 0.7110930099 D0
      W(2) = 0.2785177335 D0
      W(3) = 0.01038925650 D0

      Do i=1,3
       W(i) = fc*Z(i)/( Z(i)**4-16.0D0*Z(i)**3+72.0D0*Z(i)**2
     +                       -96.0D0*Z(i)+24.0D0 )**2
      End Do

c-------------------------
      else if(NQ.eq.4) then
c-------------------------

      Z(1) = 0.32254768961939 D0
      Z(2) = 1.74576110115835 D0
      Z(3) = 4.53662029692113 D0
      Z(4) = 9.39507091230113 D0

      W(1) = 0.6031541043 D0
      W(2) = 0.3574186924 D0
      W(3) = 0.03888790851 D0
      W(4) = 0.0005392947055 D0

      Do i=1,4
       W(i) = fc*Z(i)/(-Z(i)**5+25.0D0*Z(i)**4-200.0D0*Z(i)**3
     +               +600.0D0*Z(i)**2-600.0D0*Z(i)+120.0D0 )**2
      End Do

c--------------------
      else if(NQ.eq.5) then
c--------------------

      Z(1) =  0.26356031971814 D0
      Z(2) =  1.41340305910652 D0
      Z(3) =  3.59642577104071 D0
      Z(4) =  7.08581000585885 D0
      Z(5) = 12.64080084427580 D0

      W(1) = 0.52175 56106 D0
      W(2) = 0.39866 68111 D0
      W(3) = 0.07594 244968 D0
      W(4) = 0.00361 175868 D0
      W(5) = 0.00002 33699 7239 D0

      Do i=1,5
       W(i) = fc*Z(i)/(Z(i)**6-36.0D0*Z(i)**5+450.0D0*Z(i)**4
     +        -2400*Z(i)**3+5400.0D0*Z(i)**2-4320.0D0*Z(i)+720.0D0 )**2
      End Do

c-----------
      End if
c-----------

      Return
      End
