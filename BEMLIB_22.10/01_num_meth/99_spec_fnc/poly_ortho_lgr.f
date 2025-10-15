      function poly_ortho_lgr(n)

c----------------------------------------------
c calculate the definite integral of x^(i+j)*exp(-x)
c over the interval [0, infty)
c corresponding to the Laguerre polynomials
c----------------------------------------------

      Implicit Double Precision (a-h,o-z)

      fact=1.0D0

      if(n.gt.1) then
       Do k=1,n
        fact=fact*k
       end do
      end If

      poly_ortho_lgr =  fact

c-----
c done
c-----

      return
      end

