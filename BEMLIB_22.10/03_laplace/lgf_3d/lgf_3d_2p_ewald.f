      subroutine lgf_3d_2p_ewald
     +
     +    (a11,a12,a21,a22  ! input
     +    ,b11,b12,b21,b22
     +    ,ew,area
     +    )

c---------------------------------------------
c Compute the reciprocal lattice base vectors
c and the optimal value of the parameter xi
c---------------------------------------------

      Implicit double precision (a-h,o-z)

c----------
c constants
c----------

      pi   = 3.14159 265358 D0
      pi2  = 2.0D0*pi
      srpi = sqrt(pi)

c---------------
c unit cell area
c---------------

      area = a11*a22-a21*a12

c-----------------------------------------
c lattice base vectors in wave number space
c-----------------------------------------

      f    = pi2/area
      b11  =  f*a22
      b12  = -f*a21
      b21  = -f*a12
      b22  =  f*a11

c     write (6,*) " Reciprocal Lattice "
c     write (6,*) " ------------------ "
c     write (6,104) b11,b12
c     write (6,104) b21,b22

      ew = 0.5D0 * srpi/sqrt(area)

c     write (6,114) area,ew

c-----
c Done
c-----

 104  Format (6(1x,f10.5))
 114  Format (2x," area=",f10.5," xi_beenakker=",f10.5)

      return
      end
