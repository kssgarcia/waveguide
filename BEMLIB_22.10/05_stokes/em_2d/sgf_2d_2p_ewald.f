      subroutine sgf_2d_2p_ewald
     +
     +   (a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,area
     +   )
 
c==========================================
c BEMLIB, FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c----------------------------------------
c Computes the reciprocal lattice vectors
c          the cell area
c          xi (ew) according to Beennakker
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------
 
      pi   = 3.14159 265358 D0
      pi2  = 2.0D0*pi
      srpi = Dsqrt(pi)

c----------
c cell area
c----------

      area = Dabs(a11*a22-a12*a21)

c-----------------------------------------
c lattice base vector in wave number space 
c-----------------------------------------

      cf = pi2/area

      b11 =  cf*a22
      b12 = -cf*a21
      b21 = -cf*a12
      b22 =  cf*a11

      ew = srpi/Dsqrt(area)

c--------------------------
c     write (6,*) " Reciprocal Lattice "
c     write (6,*) " ------------------ "
c     write (6,104) b11,b12
c     write (6,104) b21,b22
c     write (6,114) area,ew
c--------------------------

c-----
c done
c-----

 104  Format (6(1x,f10.5))
 114  Format (2x," cell area = ",f10.5," ewald = ",f10.5)

      return
      end
