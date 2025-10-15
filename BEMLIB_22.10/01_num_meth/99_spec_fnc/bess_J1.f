      function bess_J1(x)

c===================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===================================

c------------------------------------------------
c This program accompanies the book:
c
c              C. Pozrikidis
c Numerical Computation in Science and Engineering
c         Oxford University Press
c------------------------------------------------

c----------------------------------------
c  Approximation of Bessel function J1
c
c  Formulas from:
c
c  Hart et al (1968) Computer Approximations,
c  John Wiley
c----------------------------------------

      Implicit Double Precision (a-h,o-z)

c---------------------------
      if(abs(x).lt.8.0) then
c---------------------------

c---
c code JONE 6040 on page 304
c formula on page 147
c---

      p00 = 72362614232.0
      p01 = -7895059235.0
      p02 =   242396853.1
      p03 =    -2972611.439
      p04 =       15704.48260
      p05 =         -30.16036606

      q00 = 144725228442.0
      q01 =   2300535178.0
      q02 =     18583304.74
      q03 =        99447.43394
      q04 =          376.9991397
      q05 =            1.0

      z = x**2

      bess_J1 = x*(p00+z*(p01+z*(p02+z*(p03+z*(p04+z*p05)))))
     +           /(q00+z*(q01+z*(q02+z*(q03+z*(q04+z*q05)))))

c---------
      else
c---------

c---
c use asymptotic series (6.8.16) on page 143
c---

      p00 =  1.0D0
      p01 =  0.00183105
      p02 = -0.00003516396496
      p03 =  0.000002457520174
      p04 = -0.000000240337019

      q00 =  0.04687499995
      q01 = -0.0002002690873
      q02 =  0.000008449199096
      q03 = -0.00000088228987
      q04 =  0.000000105787412

      t = Dabs(x)
      z = 8.0D0/t
      y = z**2

      X0 = t-2.356194491

      bess_J1 = Dsqrt(0.636619772/t)
     +     *( Dcos(X0)  *(p00+y*(p01+y*(p02+y*(p03+y*p04))))
     +       -Dsin(X0)*z*(q00+y*(q01+y*(q02+y*(q03+Y*q04))))
     +      )

c-----------
      end if
c-----------

      if(x.lt.0) bess_J1 = - bess_J1

c---
c done
c---

      return
      end
