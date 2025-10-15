      subroutine gauss_her_fun (x,y,menu)

c=========================================
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=========================================

      Implicit Double Precision (a-h,o-z)

      common/pii/pi

      If(menu.eq.1) then
       y = 1.0
      Else If(menu.eq.2) then
       y = x
      Else If(menu.eq.3) then 
       y = x**2
      Else If(menu.eq.4) then
       y = x**3
      Else If(menu.eq.5) then
       y = x**4
      Else If(menu.eq.6) then
       y = x**5
      Else If(menu.eq.7) then 
       y = x**6
      Else If(menu.eq.8) then
       y = x**7
      Else If(menu.eq.9) then
       y = x**8
      Else If(menu.eq.11) then
       y = x**10
      Else If(menu.eq.20) then
       y = 0.5D0*(x**4/16.0D0+1.0D0)*exp(x-sqrt(x**2+1.0D0))
      End If

c-----
c Done
c-----

      Return
      End
