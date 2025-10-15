      subroutine sgf_2d_2p_vvv
     +
     +   (b11,b12
     +   ,b21,b22
     +   ,max2
     +   ,ew
     +   )

c-----------------------------------------
c BEMLIB FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------------
c v_ijl is an array for summing in reciprocal space
c used to compute the stress Green's function
c
c p_i is an array for summing in reciprocal space
c used to compute the pressure Green's function
c--------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision k1,k2,kk,ks

      Dimension vxxx(-15:15,-15:15)
      Dimension vxxy(-15:15,-15:15)
      Dimension vyxy(-15:15,-15:15)

      Dimension vxyx(-15:15,-15:15)
      Dimension vxyy(-15:15,-15:15)
      Dimension vyyy(-15:15,-15:15)

      Dimension ppx (-15:15,-15:15)
      Dimension ppy (-15:15,-15:15)

      common/vvvv_2d/vxxx,vxxy,vyxy,vxyx,vxyy,vyyy,ppx,ppy

c-----------------------------------
c scan the reciprocal lattice points
c-----------------------------------

      Do i1=-max2,max2
       Do i2=-max2,max2

       k1 = i1*b11 + i2*b21
       k2 = i1*b12 + i2*b22

       ks = k1**2 + k2**2 

        if(ks.gt.0.000001) then    ! skip the zero wavenumber

        kk = Dsqrt(ks)
        t  = kk/ew
        ts = t**2

        arg = -0.25D0*ts
        EXPP = exp(arg)

        CF1 = (1.0D0+0.25D0*ts)*EXPP/kk

        k1 = k1/kk
        k2 = k2/kk

        ipress = 1
        ipress = 2

        if(ipress.eq.1) then
          CFpress = EXPP/kk
        elseif(ipress.eq.2) then
          CFpress = CF1
        end if

        ppx(i1,i2) = CFpress*k1
        ppy(i1,i2) = CFpress*k2

        vxxx(i1,i2) = (2.0*k1*k1*k1 -k1-k1)*CF1
        vxxy(i1,i2) = (2.0*k1*k1*k2 -k2   )*CF1
        vyxy(i1,i2) = (2.0*k2*k1*k2       )*CF1

        vxyx(i1,i2) = (2.0*k1*k2*k1       )*CF1
        vxyy(i1,i2) = (2.0*k1*k2*k2 -k1   )*CF1
        vyyy(i1,i2) = (2.0*k2*k2*k2 -k2-k2)*CF1

       end if

       End Do
      End Do

c-----
c done
c-----

      return
      end
