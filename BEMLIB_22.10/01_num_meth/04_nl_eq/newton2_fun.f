      subroutine newton2_fun (menu,x,f)

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c-------------------------------
c define the functions
c whose roots are to be computed
c-------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension x(2),f(2)


      if(menu==1) then
         f(1) = x(1)**2 + 2.0*x(2)**2 - 9.0
         f(2) = x(1)*x(2) -2.0
      else
         f(1) = x(1)+x(2) -3.0
         f(2) = x(1)-x(2) -2.0
      end if

c---
c done
c---

      return
      end
