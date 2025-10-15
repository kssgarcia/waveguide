      function poly_ortho_lin(n)

c----------------------------------------------
c calculates the definite integral of x^(i+j)*x
c over the interval [0,1]
c corresponding to w(x)=x
c----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      poly_ortho_lin = 1.0D0/(n+2.0D0)

c-----
c done
c-----

      return
      end
