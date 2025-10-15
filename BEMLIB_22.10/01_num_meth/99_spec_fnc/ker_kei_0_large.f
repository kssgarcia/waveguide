      subroutine ker_kei_0_large (X,KER0,KEI0)

c--------------------------------------------
c FDLIB BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c--------------------------------------------

c--------------------------------------------
c Computation of Kelvin functions ker0, kei0,
c for x>8
c--------------------------------------------

      Implicit Double Precision (a-h, o-z)

      Double Precision KER0,KEI0

      common/pip/pi,piq,pih

c----------
c constants
c----------

      two = 2.0D0
      srt = sqrt(two)

      Y  = 8.0D0/X
      Y2 = Y*Y
      Y3 = Y2*Y
      Y4 = Y2*Y2
      Y5 = Y4*Y
      Y6 = Y4*Y2

      tmp = Dsqrt(pih/X)*
     + dexp(-X/srt-0.0110486*Y+0.0000906*Y3-0.0000252*Y4
     +                       +0.0000034*Y5+0.0000006*Y6)

      AN = -X/srt-0.3926991+0.0110485*Y-0.0009765*Y2+0.0000901*Y3
     +                                 -0.0000051*Y5+0.0000010*Y6

      KER0 = tmp*cos(AN)
      KEI0 = tmp*sin(AN)

c-----
c done
c-----

      return
      end
