      subroutine lobatto (NQ,Z,W)

c======================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c======================================

c-----------------------------------------------------
c This program accompanies the book:
c             C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c         Oxford University Press
c-----------------------------------------------------

c------------------------------------------------
c Base points and weights for the Lobatto
c quadrature with NQ points
c
c Integral = sum_{i=1}^{NQ} f_i w_i

c The first base point is set at -1
c and the last base point is set at 1
c
c The intermediate NQ-2 points are located
c at the zeros of the (NQ-2)-degree Lobatto polynomial
c
c This table contains values for  NQ = 2,3,4,5,6,7,8
c
c Default value is 8
c------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension Z(20),W(20)

c-----
c trap
c-----

      if(    NQ.ne. 2
     +  .and.NQ.ne. 3
     +  .and.NQ.ne. 4
     +  .and.NQ.ne. 5
     +  .and.NQ.ne. 6
     +  .and.NQ.ne. 7
     +  .and.NQ.ne. 8
     +  ) then

        write (6,*)
        write (6,*) ' lobatto: Chosen number of base points'
        write (6,*) '          is not available; Will take 8'

        NQ = 8

      end if

c----------------------
c first and last points
c----------------------

      Z(1) = -1.0D0
      W(1)  = 2.0D0/(NQ*(NQ-1.0D0))

      Z(NQ) = 1.0D0
      W(NQ) = W(1)

c--------------------
      if(NQ.eq.3) then
c--------------------

      Z(2) = 0.0D0
      W(2) = 4.0D0/3.0D0

c-------------------------
      else if(NQ.eq.4) then
c-------------------------

      five = 5.0D0
      Z(2) = -1.0D0/Dsqrt(five)
      Z(3) = -Z(2)

      W(2) = 5.0D0/6.0D0
      W(3) = W(2)

c-------------------------
      else if(NQ.eq.5) then
c-------------------------

      Z(2) = -Dsqrt(3.0D0/7.0D0)
      Z(3) = 0.0D0
      Z(4) = -Z(2)

      W(2) = 49.0D0/90.0D0;
      W(3) = 32.0D0/45.0D0;
      W(4) = W(2);

c-------------------------
      else if(NQ.eq.6) then
c-------------------------

      Z(2) = -0.76505532392946D0
      Z(3) = -0.28523151648064D0
      Z(4) = -Z(3)
      Z(5) = -Z(2)

c     W(1) = 0.37847495629785D0
c     W(2) = 0.55485837703549D0

      Do i=2,3
      W(i) = W(1)*64.0D0/ ( (63.0D0*Z(i)**4-70.0D0*Z(i)**2+15.0D0)**2
     +      * Z(i)**2)
      End Do
      W(4) = W(3)
      W(5) = W(2)

c-------------------------
      else if(NQ.eq.7) then
c-------------------------

      Z(2) = -0.83022389627857D0
      Z(3) = -0.46884879347071D0
      Z(4) =  0.0D0
      Z(5) = -Z(3)
      Z(6) = -Z(2)

c     W(1) = 0.27682604736157
c     W(2) = 0.43174538120986
c     W(3) = 0.48761904761905

      Do i=2,4
       W(i) = W(1) * 256.0D0/( 231.0D0*Z(i)**6-315.0D0*Z(i)**4
     +             +105.0D0*Z(i)**2-5)**2
      End Do
      W(5) = W(3)
      W(6) = W(2)

c-------------------------
      else if(NQ.eq.8) then
c-------------------------

      Z(2) = -0.87174014850961D0
      Z(3) = -0.59170018143314D0
      Z(4) = -0.20929921790248D0

      Z(5) = -Z(4)
      Z(6) = -Z(3)
      Z(7) = -Z(2)

c     W(2) = 0.21070422714350
c     W(3) = 0.34112269248350
c     W(4) = 0.41245879465870

      Do i=2,4
       W(i) = W(1) * 256.0D0/ ( ( 429.0D0*Z(i)**6-693.0D0*Z(i)**4
     +              +315.0D0*Z(i)**2-35.0D0)**2 * Z(i)**2)
      End Do

      W(5) = W(4)
      W(6) = W(3)
      W(7) = W(2)

c-----------
      end if
c-----------
 
      return
      end
