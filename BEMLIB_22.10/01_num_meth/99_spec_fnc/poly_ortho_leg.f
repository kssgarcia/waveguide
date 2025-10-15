      function poly_ortho_leg(n)

c===========================================
c FDLIB
c
c calculates the definite integral of x^(i+j) 
c in the interval [-1,1]
c corresponding to the Legendre polynomials
c (unit weighting function)
c
c  prj_leg = 0          for i+j odd
c          = 2/(i+j+1)  for i+j even
c
c===========================================

      Implicit Double Precision (a-h,o-z)

      m = n/2
      index = n-2*m

      if(index.eq.0) then
        poly_ortho_leg = 2.0D0/(n+1.0D0)
      else
        poly_ortho_leg = 0.0D0
      end if

c-----
c Done
c-----

      return
      end
