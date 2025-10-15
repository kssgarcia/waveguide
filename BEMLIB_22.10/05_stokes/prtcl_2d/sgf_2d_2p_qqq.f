      subroutine sgf_2d_2p_qqq
     +
     +  (b11,b12
     +  ,b21,b22
     +  ,max2
     +  ,ew
     +  )

c-----------------------------------------
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-------------------------------------------------
c q_ij is an array for summing in reciprocal space
c
c Used to compute the velocity Green's function
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision k1,k2,kk,ks

      Dimension qxx(-15:15,-15:15)
      Dimension qxy(-15:15,-15:15)
      Dimension qyy(-15:15,-15:15)

      common/qqqq_2d/qxx,qxy,qyy

c-----------------------------------
c scan the reciprocal lattice points
c-----------------------------------

      Do i1=-max2,max2
       Do i2=-max2,max2

       k1 = i1*b11 + i2*b21
       k2 = i1*b12 + i2*b22

       ks = k1**2+k2**2

       if(ks.gt.0.000001) then   ! skip the zero wavenumber

        kk = Dsqrt(ks)
        t  = kk/ew
        ts = t*t

        arg = -0.25D0*ts
        CF = (1.0D0 + 0.25D0*ts)*exp(arg)/ks

        k1 = k1/kk
        k2 = k2/kk

        qxx(i1,i2) = (1.0D0 -k1**2)*CF
        qxy(i1,i2) =        -k1*k2 *CF
        qyy(i1,i2) = (1.0D0 -k2**2)*CF

       end if

       end Do
      end Do

c-----
c done
c-----

      return
      end
