      subroutine ell_int (RK,F,E)

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c------------------------------------------------
c This program accompanies the book:
c
c            C. Pozrikidis
c ``Numerical Computation in Science and Engineering''
c       Oxford University Press
c------------------------------------------------

c-------------------------------------------
c  Evaluation of complete elliptic integrals
c  of the first and second kind
c  using a recursion formula
c
c  SYMBOLS:
c  --------
c
c  F:	first kind
c  E:	second kind
c  acc: specified accuracy
c-------------------------------------------

      Implicit Double Precision (a-h,o-z)

      Parameter (acc=0.000000001)

c----------
c constants
c----------

      pi  = 3.14159 265358  D0
      pih = 0.5D0*pi
      
c----------
c launching
c----------

      F  = pih
      E  = 1.0D0
      G  = 1.0D0
      B  = RK
 
c     Do while (abs(D).GT.acc)

  1   Continue

      C = Dsqrt(1.0D0-B*B)
      B = (1.0D0-C)/(1.0D0+C)
      D = F*B
      F = F+D
      G = 0.5D0*G*B
      E = E+G

      If(abs(D).gt.acc) Go to 1

c     End Do

      E = F*(1.0D0-0.5D0*RK*RK*E)

c-----
c Done
c-----
 
      Return
      End
