      subroutine sgf_3d_3p_ewald
     +
     +  (a11,a12,a13
     +  ,a21,a22,a23
     +  ,a31,a32,a33
     +  ,b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,ew
     +  ,vlm
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c-------------------------------------
c computes the reciprocal base vectors
c and the parameter xi according to Beenaker
c-------------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi   = 3.14159 265358 D0
      pi2  = 2.0D0*pi
      srpi = Dsqrt(pi)

      oot  = 1.0D0/3.0D0

c------------
c cell volume
c------------

      vlm = a31*(a12*a23-a22*a13)
     +    + a32*(a13*a21-a23*a11)
     +    + a33*(a11*a22-a21*a12)

c-----------------------------------------
c lattice base vectors in wave number space
c-----------------------------------------

      fc = pi2/vlm

      b11 = fc*(a22*a33-a23*a32)
      b12 = fc*(a23*a31-a21*a33)
      b13 = fc*(a21*a32-a22*a31)

      b21 = fc*(a32*a13-a33*a12)
      b22 = fc*(a33*a11-a31*a13)
      b23 = fc*(a31*a12-a32*a11)

      b31 = fc*(a12*a23-a13*a22)
      b32 = fc*(a13*a21-a11*a23)
      b33 = fc*(a11*a22-a12*a21)

      ew = srpi/vlm**oot

c-----------------------------------------
c     write (6,*) " Reciprocal Lattice "
c     write (6,*) " ------------------ "
c     write (6,104) b11,b12,b13
c     write (6,104) b21,b22,b23
c     write (6,104) b31,b32,b33
c     write (6,114) vlm,ew
c-----------------------------------------

c-----
c done
c-----

 104  Format (6(1x,f10.5))
 114  Format (2x," vlm=",f10.5," xi_beenakker=",f10.5)

      return
      end
