      subroutine ell_int (RK2,F,E)

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-------------------------------------------
c  EVALUATION OF COMPLETE ELLIPTIC INTEGRALS
c  OF THE FIRST AND SECOND KIND
c
c  F:   first kind
c  E:   second kind
c-------------------------------------------

      Implicit Double Precision (A-H,O-Z)

      parameter (acc=0.000000001)

c----------
c constants
c----------

      pi  = 3.14159 265358 D0
      pih = 0.5D0*pi

c---------------------
c compute by iteration
c---------------------

      RK = sqrt(RK2)
      F  = pih
      E  = 1.0D0
      G  = 1.0D0
      B  = RK

  1   Continue

c     Do 99 while (abs(D).GT.acc)

      C = sqrt(1.0D0-B**2)
      B = (1.0D0-C)/(1.0D0+C)
      D = F*B
      F = F+D
      G = 0.5D0*G*B
      E = E+G

c     End Do

      If(abs(D).gt.acc) Go to 1

      E=F*(1.0D0-0.5D0*RK2*E)

c-----
c Done
c-----

      return
      end
