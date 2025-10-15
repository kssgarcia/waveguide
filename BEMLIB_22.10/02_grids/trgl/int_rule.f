      subroutine int_rule (mint,Nbase,xi,eta,w)

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c---------------------------------
c Generate base points and weights
c for a triangle integration rule
c--------------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  xi(250000)
      Dimension eta(250000)
      Dimension   w(250000)    ! integration weights

c--------
c base points
c--------

      Dxi = 1.0D0/mint

       Ic =0

       Do i=1,mint+1    ! integation rule
        Do j=1,mint+2-i   ! integation rule
         Ic=Ic+1
         xi(Ic)  = (i-1.0D0)*Dxi  ! integation rule
         eta(Ic) = (j-1.0D0)*Dxi  ! integation rule
        End Do
      End Do

      Nbase = Ic

c--------
c prepare
c--------

      fc = Dxi**2

      wsum = 0.0D0
      Ic=0

c-----
c scan
c-----

      Ic = Ic+1
      w(Ic) = fc/4.0D0
      wsum = wsum+w(Ic)

      Do j=2,mint-1
        Ic = Ic+1
        w(Ic) = 0.5D0*fc
        wsum = wsum+w(Ic)
      End Do

      Ic = Ic+1
      w(Ic) = 5.0D0/12.0D0 * fc
      wsum = wsum+w(Ic)

      Ic = Ic+1
      w(Ic) = 1.0D0/6.0D0 * fc
      wsum = wsum+w(Ic)

      Do i=2,mint-1

       Ic = Ic+1
       w(Ic) = 0.5D0 * fc
       wsum = wsum+w(Ic)

       Do j=2,mint-i
        Ic = Ic+1
        w(Ic) = fc
        wsum = wsum+w(Ic)
       End Do

       Ic = Ic+1
       w(Ic) = 11.0D0/12.0D0 * fc
       wsum = wsum+w(Ic)
       Ic = Ic+1
       w(Ic) = 7.0D0/12.0D0 * fc
       wsum = wsum+w(Ic)

      End Do
c---

      Ic = Ic+1
      w(Ic) = 5.0D0/12.0D0 * fc
      wsum = wsum+w(Ic)

      Ic = Ic+1
      w(Ic) = 7.0D0/12.0D0 * fc
      wsum = wsum+w(Ic)

c---

      Ic = Ic+1
      w(Ic) = 1.0D0/6.0D0 * fc
      wsum = wsum+w(Ic)

c-----
c done
c-----

c     write (6,*) Ic,(mint+1)*(mint+2)/2,wsum

      return
      end
