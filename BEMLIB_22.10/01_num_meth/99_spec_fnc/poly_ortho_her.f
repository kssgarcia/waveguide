      function poly_ortho_her(n)

c=========================================
c calculates the definite integral of x^n*exp(-x^2)
c over the interval [0, infty]
c corresponding to the Hermite polynomials
c=========================================

      Implicit Double Precision (a-h,o-z)

      common/pii/pi,srpi

      m = n/2
      index = n-2*m

      if(index.eq.0) then
         fact = 1.0
         Do i=1,m
          fact=(2*i-1)*fact
         End Do
         poly_ortho_her = fact*srpi/2**m
      else
         poly_ortho_her = 0.0D0
      end if

c-----
c done
c-----

      return
      end
