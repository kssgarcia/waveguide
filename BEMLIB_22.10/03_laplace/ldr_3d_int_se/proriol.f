      subroutine proriol(k,l,xi,eta,pror)

c========================================
c FDLIB, BEMLIB
c
c evaluation of the kl Proriol polynomial
c========================================

      Implicit Double Precision (a-h,o-z)

      xip = 2*xi/(1-eta)-1
      etap = 2*eta-1
      pror = jacobi(0,0,k,xip) * (1-eta)**k
     +     * jacobi(2*k+1,0,l,etap)

c-----
c done
c-----

      return
      end
