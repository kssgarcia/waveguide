      function poly_ortho_rad(n)

c------------------------------------------------------
c calculate the definite integral of x^(i+j) * (1.0+x)
c over the interval [-1,1]
c corresponding to the Radau polynomials
c------------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      m = n/2
      index = n-2*m

      if(index.eq.0) then
         poly_ortho_rad = 2.0D0/(n+1.0D0)
      else
         poly_ortho_rad = 2.0D0/(n+2.0D0)
      end if

c-----
c Done
c-----

      return
      end

