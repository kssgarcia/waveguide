      subroutine lobatto_dr_fun (x,f,menu)

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

      Implicit Double Precision (a-h,o-z)

      common/pii/pi

      if(menu.eq.1) then
       f = 1.0D0
      else if(menu.eq.2) then
       f = Dexp(pi*x/2.0D0)
      else if(menu.eq.3) then
       f = x*x*x*x
      end if

c-----
c done
c-----

      return
      end
