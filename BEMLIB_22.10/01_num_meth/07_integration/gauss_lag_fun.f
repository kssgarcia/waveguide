      subroutine gauss_lag_fun (x,y,menu)

c==========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c--------------------------------------------
c Function definition for the Gauss--Laguerre
c quadrature
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)

      common/pii/pi

      if(menu.eq.1) then
       y = 1.0
      else if(menu.eq.2) then
       y = x
      else if(menu.eq.3) then 
       y = x**2
      else if(menu.eq.4) then
       y = x**3
      else if(menu.eq.5) then
       y = x**4
      else if(menu.eq.6) then
       y = x**5
      else if(menu.eq.7) then 
       y = x**6
      else if(menu.eq.8) then
       y = x**7
      else if(menu.eq.9) then
       y = x**8
      else if(menu.eq.10) then
       y = x**9
      else if(menu.eq.11) then
       y = x**10
      else if(menu.eq.20) then
       y = 0.5*(x**4/16.0+1.0)*exp(x-sqrt(x**2+1))
      else if(menu.eq.25) then
       y = 0.5*(x**4/16.0+1.0)*exp(-10/x**2)
      end if

c-----
c Done
c-----

      Return
      End
