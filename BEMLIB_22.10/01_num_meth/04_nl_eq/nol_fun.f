      subroutine nol_fun (menu,x,f)

c=======================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c=======================================

c--------------------------------
c Define the functions comprising
c the nonlinear system
c
c menu = 1, 2:
c
c   x(1) = xi,  x(2) = eta
c--------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension x(20),f(20)

      common/pii/pi
      common/menu12/al
      common/menu5/a,b,c

c-----------------------
      if(menu.eq.1) then   ! anti-symmetric Stokes flow
c-----------------------

        q     = 2.0D0*al
        RK    = sin(q)/q
        coshy = 0.5D0*(dexp(x(2))+dexp(-x(2)))
        sinhy = 0.5D0*(dexp(x(2))-dexp(-x(2)))
        f(1)  = dsin(x(1))*coshy + RK*x(1)
        f(2)  = dcos(x(1))*sinhy + RK*x(2)

c-----------------------
      elseif(menu.eq.2) then  ! symmetric Stokes flow
c-----------------------

        q     = 2.0D0*al
        RK    = sin(q)/q
        coshy = 0.5D0*(dexp(x(2))+dexp(-x(2)))
        sinhy = 0.5D0*(dexp(x(2))-dexp(-x(2)))
        f(1)  = dsin(x(1))*coshy - RK*x(1)
        f(2)  = dcos(x(1))*sinhy - RK*x(2)

c-----------------------
      elseif(menu.eq.3) then
c-----------------------

        f(1) = (2.0*x(1)-x(2))**2
        f(2) = (x(1)-1.0)*(x(2)-2.0)

c-----------------------
      else if(menu.eq.4) then
c-----------------------

        f(1) = 3.0*x(1)-cos(x(2)*x(3))-0.5
        f(2) = x(1)**2-81.0*(x(2)+0.1)**2+sin(x(3))+1.06
        f(3) = exp(-1.0*x(1)*x(2))+20.0*x(3)+(10.0*pi-3.0)/3.0

c-----------------------
      elseif(menu.eq.5) then
c-----------------------

        ab = a+b
        ac = a+c
        bc = b+c
        f(1) = ab-ac*cos(x(1))-bc*cos(x(2))
        f(2) =    ac*sin(x(1))-bc*sin(x(2))

c-----------------------
      elseif(menu.eq.6) then
c-----------------------

        f(1) = (0.5*x(1)-x(2))*(0.5*x(1)+x(2))
     +       - 0.312 * (1.0-x(1)-x(2))**2
        f(2) = x(2)*(0.5*x(1)+x(2))
     +       - 0.480 * (1.0-x(1)-x(2))*(0.50*x(1)-x(2))

c-----------------------
      elseif(menu.eq.7) then
c-----------------------

        f(1) = (x(1)-1.0)**2+(x(2)-1.0)**2-2.0
        f(2) = x(1)*x(2)-4.0

c-----------------------
      elseif(menu.eq.8) then
c-----------------------

        f(1) = (x(1)-2.0)**2+(x(2)-3.0)**2+(x(1)-2.1)*(x(2)-3.1)-2.81
        f(2) = 10.0*exp(-x(1)) + 5.0*exp(1-x(2))-0.7468

c-----------
      end if
c-----------

      return
      end
