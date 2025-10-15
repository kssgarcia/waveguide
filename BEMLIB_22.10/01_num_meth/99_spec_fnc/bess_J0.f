      function bess_J0(x)

c=========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

c------------------------------------------------
c This program accompanies the book:
c        C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c      Oxford University Press
c------------------------------------------------

c-----------------------------------------
c  Approximation of Bessel function J0
c
c  Formulas from:
c
c  Hart et al, (1968)
c  Computer Approximations, John Wiley
c-----------------------------------------

      Implicit Double Precision (a-h,o-z)

c---------------------------
      if(abs(x).lt.8.0) then
c---------------------------

c---
c code JZERO 5838 on page 298
c formula on page 147
c---

      p00 =  57568490574.0
      p01 = -13362590354.0
      p02 =    651619640.7
      p03 =    -11214424.18
      p04 =        77392.33017
      p05 =         -184.9052456

      q00 =  57568490411.0
      q01 =   1029532985.0
      q02 =      9494680.718
      q03 =        59272.64853
      q04 =          267.8532712
      q05 =            1.0D0

      z = x**2

      bess_J0 = (p00+z*(p01+z*(p02+z*(p03+z*(p04+z*p05)))))
     +         /(q00+z*(q01+z*(q02+z*(q03+z*(q04+z*q05)))))

c---------
      else
c---------

c---
c use asymptotic series (6.8.16) on page 143
c---

      p00 =  1.0D0
      p01 = -0.001098628627
      p02 =  0.00002734510407
      p03 = -0.000002073370639
      p04 =  0.0000002093887211

      q00 = -0.01562499995
      q01 =  0.0001430488765
      q02 = -0.000006911147651
      q03 =  0.0000007621095161
      q04 = -0.0000000934945152

      t  = abs(x)
      z  = 8.0D0/t
      y  = z**2
      X0 = t-0.785398164

      bess_J0 
     +    = sqrt(0.636619772/t)
     +   *( Dcos(X0)*  (p00+y*(p01+y*(p02+y*(p03+y*p04))))
     +     -Dsin(X0)*z*(q00+y*(q01+y*(q02+y*(q03+y*q04))))
     +    )

c-----------
      end if
c-----------

c---
c done
c---

      return
      end
