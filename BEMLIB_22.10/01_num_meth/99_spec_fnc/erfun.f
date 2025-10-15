      function erfun(x)

c------------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c------------------------------------------

c----------------
c error function 
c----------------

      Implicit Double Precision (a-h,o-z)

      t = 1.0D0/(1.0D0+0.3275911*Dabs(x))

c------------------------------
c  complementary error function 
c------------------------------

      erfunc = exp(-x*x)*t*( 0.254829592
     +                  +t*(-0.284496736
     +                  +t*( 1.421413741
     +                  +t*(-1.453152027
     +                  +t*  1.061405429))))

      if(x.lt.0.0) erfunc = 2.0D0-erfunc

c----------------
c  error function 
c----------------
 
      erfun = 1.0D0-erfunc

c-----
c done
c-----

      return
      end
