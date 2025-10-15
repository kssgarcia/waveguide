      subroutine trgl_rule 
     +
     +  (mint
     +  ,nbase
     +  ,xi,eta,w
     +  )

c======================================
c base points and weights
c for triangle integration rule
c in xi and eta
c
c     |\
c eta | \
c     |  \
c     -----
c      xi
c
c mint:  descretization level
c nbase: number of base points
c======================================

      Implicit Double Precision (a-h,o-z)

      Dimension  xi(250000)
      Dimension eta(250000)
      Dimension   w(250000)    ! integration weights

c------------
c base points
c------------

      Dxi = 1.0D0/mint

       Ic = 0

       Do i=1,mint+1  
        Do j=1,mint+2-i
         Ic=Ic+1
         xi(Ic)  = (i-1.0D0)*Dxi
         eta(Ic) = (j-1.0D0)*Dxi
        End Do
      End Do

      nbase = Ic

c--------
c prepare
c--------

      fc = Dxi*Dxi

      wsum = 0.0D0

c-----
c scan
c-----

      Ic = 1
      w(Ic)=fc/4.0D0
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
