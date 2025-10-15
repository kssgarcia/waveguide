      subroutine newton1_fun (menu,x,f)

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c---------------------------------------------------
c define the function whose roots are to be computed
c---------------------------------------------------

      Implicit Double Precision (a-h,o-z)

      common/menu/a,b,OH,OOT,RT,P,T,al,eps_con,rl,rlp
      common/pii/pi
 
c------------------------
      if(menu.eq.1) then
c------------------------

      R = sqrt(1.0-1.0/(x**2))
      f = A- sqrt(1.0+x*asin(R)/R) / (RT*X**OOT)

c------------------------
      elseif(menu.eq.2) then
c------------------------

      X6 = 30/(3.1-X)
      X7 = (45-2.361378*X6)/(1-X)
      X8 = (25-0.165610*X6)/(1-X)
      f = 33*X6+33*X7-67*X8

c------------------------
      elseif(menu.eq.3) then
c------------------------

      alpha = 2.0
      beta  = 4.5
      ca    = x
      cb    = alpha-2.0*(1.0-ca)
      f     = 1.0-ca-beta*ca*cb**2

c------------------------
      elseif(menu.eq.4) then
c------------------------

      R = 82.06
      a = 1418.612
      b = 18.4

      f = 0.000001*P*(x**2-b**2)*x-R*T*x*(x+b)*0.000001
     +   +a*(X-b)/sqrt(T)

c------------------------
      elseif(menu.eq.5) then
c------------------------

      f = exp(x)+cos(x)-2.52

c------------------------
      elseif(menu.eq.6) then
c------------------------

      f = 3*x**2 + eps_con*5.0*log(abs(pi-x))+1.0
     +   +3.0*(eps_con-1.0)*x**3

c------------------------
      elseif(menu.eq.7) then
c------------------------

      f = 9.0-2.0*x**2+5.0*x**3-x**4

c------------------------
      elseif(menu.eq.8) then
c------------------------

      f = log(x)+a*x-b

c-----------------------------
      elseif(menu.eq.9) then
c-----------------------------

      gm = x
      fc = rl/(2.0D0-rl)
      bt = atan(fc*tan(gm))
      f =    sin(  rl  *al-bt)/cos(bt)
     +   -fc*sin((rl-2)*al-gm)/cos(gm)

c-----------------------------
      elseif(menu.eq.10) then
c-----------------------------

      f = log(x)+3-3.1*x*x

c-----------------------------
      elseif(menu.eq.11) then
c-----------------------------

      f = x*log(abs(x))

c-----------------------------
      elseif(menu.eq.101) then
c-----------------------------

      f = sin(2.0D0*al*(x-1.0D0))-(1.0D0-x)*sin(2.0D0*al)

c-----------------------------
      elseif(menu.eq.102) then
c-----------------------------

      f = sin(2.0*al*(x-1.0D0))-(x-1.0D0)*sin(2.0D0*al)
      f = f/(x-2.0D0)

c-----------------------------
      elseif(menu.eq.110) then
c-----------------------------

      rj0 = bessel_J0(x)
      rj1 = bessel_J1(x)
      f = 1.0D0 + (rj1/rj0-2.0D0/x)*rj1/rj0
      f = f-2.0D0*rlp

c-----------
      end if
c-----------

c---
c done
c---

      return
      end
