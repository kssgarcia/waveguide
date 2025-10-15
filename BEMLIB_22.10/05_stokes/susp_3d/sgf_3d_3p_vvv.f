      subroutine sgf_3d_3p_vvv
     +
     +   (b11,b12,b13
     +   ,b21,b22,b23
     +   ,b31,b32,b33
     +   ,Max2
     +   ,ew
     +   )

c============================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licencing agreement
c============================================

c--------------------------------------------
c v(i,j,k):
c
c array for summing in reciprocal space
c used to compute the stress Green's function
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Double Precision k1,k2,k3,kk,ks

      Dimension vxxx(-9:9,-9:9,-9:9)
      Dimension vxxy(-9:9,-9:9,-9:9)
      Dimension vxxz(-9:9,-9:9,-9:9)
      Dimension vyxy(-9:9,-9:9,-9:9)
      Dimension vyxz(-9:9,-9:9,-9:9)
      Dimension vzxz(-9:9,-9:9,-9:9)

      Dimension vxyx(-9:9,-9:9,-9:9)
      Dimension vxyy(-9:9,-9:9,-9:9)
      Dimension vxyz(-9:9,-9:9,-9:9)
      Dimension vyyy(-9:9,-9:9,-9:9)
      Dimension vyyz(-9:9,-9:9,-9:9)
      Dimension vzyz(-9:9,-9:9,-9:9)

      Dimension vxzx(-9:9,-9:9,-9:9)
      Dimension vxzy(-9:9,-9:9,-9:9)
      Dimension vxzz(-9:9,-9:9,-9:9)
      Dimension vyzy(-9:9,-9:9,-9:9)
      Dimension vyzz(-9:9,-9:9,-9:9)
      Dimension vzzz(-9:9,-9:9,-9:9)

      Dimension ppx(-9:9,-9:9,-9:9)
      Dimension ppy(-9:9,-9:9,-9:9)
      Dimension ppz(-9:9,-9:9,-9:9)

c-----------------
c common statement
c-----------------

      common/vvvv_3d/vxxx,vxxy,vxxz,vyxy,vyxz,vzxz
     +              ,vxyx,vxyy,vxyz,vyyy,vyyz,vzyz
     +              ,vxzx,vxzy,vxzz,vyzy,vyzz,vzzz
     +              ,ppx,ppy,ppz

c-----------------------------------
c scan the wave number lattice points
c------------------------------------

      Do i1 = -Max2,Max2
      Do i2 = -Max2,Max2
      Do i3 = -Max2,Max2

       k1 = i1*b11 + i2*b21 + i3*b31
       k2 = i1*b12 + i2*b22 + i3*b32
       k3 = i1*b13 + i2*b23 + i3*b33

       ks = k1**2 + k2**2 +k3**2

c---
       if(ks.gt.0.000001) then   ! skip the zero wavenumber
c---

       kk = Dsqrt(ks)
       t  = kk/ew
       t2 = t**2
       t4 = t2**2

       arg = -0.25D0*t2
       EXPP = Dexp(arg)

       k1 = k1/kk    ! normalize
       k2 = k2/kk
       k3 = k3/kk

c---
c pressure matrix
c---

       CF2 = EXPP/kk
       CF1 = 2.0D0 * (1.0D0 + 0.25D0*t2 + 0.125D0*t4) * CF2

c---
       Ipress=3
       Ipress=2

       if(Ipress.eq.1) then
         CFP =  (1.0D0 - t2/6.0D0 - t4/12.0D0) * CF2
       else if(Ipress.eq.2) then
         CFP =  (1.0D0 + 0.25D0* t2 + 0.125D0*t4) * CF2
       else if(Ipress.eq.3) then
         CFP = CF2
       end if
c---

       ppx(i1,i2,i3) = CFP*k1
       ppy(i1,i2,i3) = CFP*k2
       ppz(i1,i2,i3) = CFP*k3

       vxxx(i1,i2,i3) = CF1*k1*k1*k1 - CF2*(k1+k1+k1)
       vxxy(i1,i2,i3) = CF1*k1*k1*k2 - CF2* k2
       vxxz(i1,i2,i3) = CF1*k1*k1*k3 - CF2* k3
       vyxy(i1,i2,i3) = CF1*k2*k1*k2 - CF2* k1
       vyxz(i1,i2,i3) = CF1*k2*k1*k3 
       vzxz(i1,i2,i3) = CF1*k3*k1*k3 - CF2* k1

       vxyx(i1,i2,i3) = CF1*k1*k2*k1 - CF2* k2
       vxyy(i1,i2,i3) = CF1*k1*k2*k2 - CF2* k1
       vxyz(i1,i2,i3) = CF1*k1*k2*k3
       vyyy(i1,i2,i3) = CF1*k2*k2*k2 - CF2*(k2+k2+k2)
       vyyz(i1,i2,i3) = CF1*k2*k2*k3 - CF2* k3
       vzyz(i1,i2,i3) = CF1*k3*k2*k3 - CF2* k2

       vxzx(i1,i2,i3) = CF1*k1*k3*k1 - CF2* k3
       vxzy(i1,i2,i3) = CF1*k1*k3*k2
       vxzz(i1,i2,i3) = CF1*k1*k3*k3 - CF2* k1
       vyzy(i1,i2,i3) = CF1*k2*k3*k2 - CF2* k3
       vyzz(i1,i2,i3) = CF1*k2*k3*k3 - CF2* k2
       vzzz(i1,i2,i3) = CF1*k3*k3*k3 - CF2*(k3+k3+k3)

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
