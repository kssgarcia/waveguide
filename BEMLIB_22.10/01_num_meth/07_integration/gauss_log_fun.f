      subroutine gauss_log_fun (x,y,menu)

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

      Implicit Double Precision (a-h,o-z)

      common/pii/pi

      if(menu.eq.1) then
       y = 1.0D0
      else if(menu.eq.2) then
       y = -Dexp(pi*x/2.0D0)
      end if

c----
c Done
c----

      Return
      End
