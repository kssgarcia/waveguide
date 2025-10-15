      subroutine sgf_2d_1p_ww_fga (H,om,sh,www)

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

      Implicit Double Precision (a-h,o-z)

      tmp = 4.0D0*H-sh
      www = (sh*(exp(-sh*om) - exp((sh-8.0D0*H)*om))
     +       + (sh-4.0D0*H)*(exp((sh-4.0D0*H)*om)
     +      - exp(-(sh+4.0D0*H)*om)))
     +      /(1.0D0-exp(-4.0D0*H*om))**2

c-----
c done
c-----

      return
      end
