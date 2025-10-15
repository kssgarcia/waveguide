      subroutine sgf_3d_2p_ewald
     +
     +   (a11,a12,a21,a22
     +   ,b11,b12,b21,b22
     +   ,ew,area
     +   )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c Computes the reciprocal lattice base vectors
c and the parameter xi according to Beenakker
c
c SYMBOLS:
c --------
c
c ew: ewald splitting constant xi
c-----------------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi   = 3.14159 265358 D0
      pi2  = 2.0D0*pi
      srpi = sqrt(pi)

c----------
c cell area
c----------

      area = a11*a22-a21*a12

c------------------------------------------
c lattice base vectors in wave number space
c------------------------------------------

      fc = pi2/area

      b11 =  fc*a22
      b12 = -fc*a21
      b21 = -fc*a12
      b22 =  fc*a11

      ew = srpi/Dsqrt(area)

c--------------------------
c     write (6,*) " Reciprocal Lattice "
c     write (6,*) " ------------------ "
c     write (6,104) b11,b12
c     write (6,104) b21,b22
c     write (6,114) area,ew
c--------------------------

 104  Format (6(1x,f10.5))
 114  Format (2x," area=",f10.5," recommended xi=",f10.5)

      return
      end
