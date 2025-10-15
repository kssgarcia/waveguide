      subroutine velocity (menu,x,t,f)

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c---------------------------------------
c Evaluation of the phase-space velocity
c---------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension x(10),f(10)

c-----
c common blocks
c-----

      common/menu1/a
      common/menu2/rk,r,b
      common/menu3/rnu
      common/menu4/rlam,pinf,epsilon

      common/pii/pih,pi

c-----------------------
      If(menu.eq.1) then
c-----------------------

       f(1) = a*x(1)

c----------------------------
      Else If(menu.eq.2) then
c----------------------------

       f(1) = - rk*(x(1)-x(2))
       f(2) = - x(1)*x(3)+r*x(1)-x(2)
       f(3) =   x(1)*x(2)-b*x(3)

c----------------------------
      Else If(menu.eq.3) then
c----------------------------

       f(1) = 1.0-0.4*t**rnu * (1.0-t)**rnu

c----------------------------
      Else If(menu.eq.4) then
c----------------------------

      rlam3 = 3.0D0*rlam

      x1s = x(1)**2
      x1q = x(1) * x(1)**rlam3
      x2s = x(2)**2

      f(1) = x(2)

      f(2) = -1.5D0*x2s/x(1) 
     +       - 2.0D0/x1s
     +       + 2.0D0/x1q
     +       + (1.0D0/x1q-1.0D0/x(1))*pinf
     +       -  pinf*epsilon/x(1)

c----------------------------
      Else If(menu.eq.5) then
c----------------------------

       f(1) = x(2)
       f(2) = x(1)*x(2) + Dsin(pih*x(1))

c----------------------------
      Else If(menu.eq.6) then
c----------------------------

       rs = x(1)**2+x(2)**2+x(3)**2
       r  = Dsqrt(rs)
       rc = r*rs 
       rq = rs*rs 
       rp = rc*rs 
       rx = rc*rc
       re = rp*rc
       rn = re*r
       rt = re*rs
       rl = rt*r
       rv = rt*rs

       If(r.ge.2.5D0) then
        A = 5.0D0/rc-8.0D0/rp+25.0D0/rx-35.0D0/re+125.0D0/rn
     +       -102.0D0/rt+12.5D0/rl+430.0D0/rv
        B = (16.0D0/rp+10.0D0/re-36.0D0/rt-25.0D0/rl-36.0D0/rv)/3.0D0
       Else If(r.lt.2.5D0.and.r.gt.2.01D0) then
        A = -4.3833D0+17.7176D0/r+14.8204D0/rs-92.4471D0/rc
     +      -46.3151D0/rq+232.2304D0/rp
        B = - 3.1918D0+12.3641D0/r+11.4615D0/rs-65.2926/rc
     +      -36.4909D0/rq+154.8074D0/rp
       Else
        C = - Dlog(r-2.0D0)
        Cs = C**2
        A = 16.3096D0/r-7.1548D0
        B = 2.0D0*(0.4056D0*Cs+1.49681D0*C-1.9108D0)/
     +     (r*(Cs+6.04250D0*C+6.32549D0))
       End If

       e = x(1)*x(2)*(B-A)/rs

       f(1) = x(2) +e*x(1)-0.50D0*B*x(2)
       f(2) =       e*x(2)-0.50D0*B*x(1)
       f(3) =       e*x(3)

c-----------
      End If
c-----------

c-----
c Done
c-----

      Return
      End
