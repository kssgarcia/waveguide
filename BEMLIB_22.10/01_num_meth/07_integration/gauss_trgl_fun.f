      subroutine gauss_trgl_fun (xi,eta,f,menu)

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

      Implicit Double Precision (a-h,o-z)

      common/pii/pi

      if(menu.eq.1) then
       f = 1.0D0
      else if(menu.eq.2) then
       f = xi
      else if(menu.eq.3) then
       f = eta
      else if(menu.eq.4) then
       f = xi**2
      else if(menu.eq.5) then
       f = eta**2
      else if(menu.eq.6) then
       f = xi*eta
      else if(menu.eq.7) then
       f = xi**3
      else if(menu.eq.8) then
       f = xi**2*eta
      end if

c-----
c Done
c-----

      Return
      End
