      function poly_ortho_log(n)

c-------------------------------------------------
c calculate the definite integral of x^n*ln(x)
c over the interval [0,1]
c-------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      poly_ortho_log = 1.0D0/(n+1.0D0)**2

c-----
c Done
c-----

      return
      end
