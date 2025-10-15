      subroutine sgf_3d_3p_qqq
     +
     +  (b11,b12,b13
     +  ,b21,b22,b23
     +  ,b31,b32,b33
     +  ,Max2
     +  ,ew
     +  )

c-----------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c-----------------------------------------
c q(i,j): 
c
c array for summing in reciprocal space
c used to compute the velocity
c Green's function
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision k1,k2,k3,kk,ks

      Dimension qxx(-9:9,-9:9,-9:9)
      Dimension qxy(-9:9,-9:9,-9:9)
      Dimension qxz(-9:9,-9:9,-9:9)
      Dimension qyy(-9:9,-9:9,-9:9)
      Dimension qyz(-9:9,-9:9,-9:9)
      Dimension qzz(-9:9,-9:9,-9:9)

      common/qqqq_3d/qxx,qxy,qxz,qyy,qyz,qzz

c-----------------------------------
c scan the reciprocal lattice points
c-----------------------------------

      Do i1 = -Max2,Max2
      Do i2 = -Max2,Max2
      Do i3 = -Max2,Max2

       k1 = i1*b11 + i2*b21 + i3*b31
       k2 = i1*b12 + i2*b22 + i3*b32
       k3 = i1*b13 + i2*b23 + i3*b33

       ks = k1**2 + k2**2 + k3**2

c---
       if(ks.gt.0.00000001) then  ! skip the zero wavenumber
c---
 
        kk = sqrt(ks)
        t  = kk/ew
        t2 = t*t
        t4 = t2*t2

        arg = -0.25D0 * t2
        EXPP = dexp(arg)
        CF = (1.0D0 + 0.25D0*t2 + 0.125D0*t4)*EXPP/ks

        k1 = k1/kk
        k2 = k2/kk
        k3 = k3/kk

        qxx(i1,i2,i3) = (1.0D0 -k1*k1)*CF
        qxy(i1,i2,i3) =        -k1*k2 *CF
        qxz(i1,i2,i3) =        -k1*k3 *CF
        qyy(i1,i2,i3) = (1.0D0 -k2*k2)*CF
        qyz(i1,i2,i3) =        -k2*k3 *CF
        qzz(i1,i2,i3) = (1.0D0 -k3*k3)*CF

c---
      end if
c---
      
c---
      end do
      end do
      end do
c---

c-----
c done
c-----

      return
      end
