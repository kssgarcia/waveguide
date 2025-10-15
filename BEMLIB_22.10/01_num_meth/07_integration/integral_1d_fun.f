      subroutine integral_1d_fun (x,y,menu)

c===========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c===========================================

c-------------------------------
c compute (define) the integrand
c-------------------------------

      Implicit Double Precision (a-h,o-z)

      common/ccc/rks,cnt
      common/vari/m
      common/varr/c
      common/pii/pi
      common/option72/alpha72,p72,q72
      common/option73/alpha73,pmq73

c-----------------------
      if(menu.eq.1) then
c-----------------------

        if(abs(x).lt.0.00001) then
          y = 0
        else
          xs = x**2
          y = -2.0D0*xs*(log(x)-1.0)*exp(-xs)
        end if

c-----------------------
      elseif(menu.eq.2) then
c-----------------------

        y = 3.0

c-----------------------
      elseif(menu.eq.3) then
c-----------------------

        y = Dlog(1.0-0.99*Dcos(x*pi/2.0)**2)

c-----------------------
      elseif(menu.eq.4) then
c-----------------------

        y = exp(x)

c-----------------------
      elseif(menu.eq.5) then
c-----------------------

        y = x**3
        y = Dsqrt(y)

c-----------------------
      elseif(menu.eq.6) then
c-----------------------

        y = exp(-0.25*(1.0+x)**2 )

c-----------------------
      elseif(menu.eq.7) then
c-----------------------

        if(x<0.00001) then
c         y = 0
          y = 1.0/4.0
        else
c         y = (1.0D0-exp(-x))/sqrt(x**3)
          y = (1.0D0-exp(-x*x/4.0D0))/x**2
        end if

c-----------------------
      elseif(menu.eq.20) then
c-----------------------
      
        y = sqrt(1.0-rks*sin(x)**2)

c-----------------------
      elseif(menu.eq.21) then
c-----------------------

        y = log(abs(cnt-sin(x)))

c-----------------------
      elseif(menu.eq.22) then
c-----------------------

        y = 2.0*(1.0-exp(-x*x/4))/(x*x)

c-----------------------
      elseif(menu.eq.30) then
c-----------------------

        y = (1.0-cos(m*x))/sin(0.5*x)**2

c-----------------------
      elseif(menu.eq.31) then
c-----------------------

        y = 1.0/(m**2+sin(0.5*x)**2)

c-----------------------
      elseif(menu.eq.32) then
c-----------------------

        y = x**2/sqrt(2.0-x**2)

c-----------------------
      elseif(menu.eq.41) then
c-----------------------

        y = (cos(2.0D0*m*x)-cos(2.0D0*(m-1.0D0)*x))/sin(x)

c-----------------------
      elseif(menu.eq.42) then
c-----------------------

        sn = sin(x)
        den = sn*sqrt(1+sn**2)
        y = (cos(2.0D0*m*x)-cos(2.0D0*(m-1.0D0)*x))/den

c-----------------------
      elseif(menu.eq.45) then
c-----------------------

        cs = cos(2.0D0*x)
        den = sqrt(cs**2-8.0D0*cs+7.0D0)
        y = (1.0D0-cos(2.0D0*m*x))/den

c-----------------------
      elseif(menu.eq.60) then
c-----------------------

        y = cos(c*x)*(1.0D0-x**2)

c-----------------------
      elseif(menu.eq.61) then
c-----------------------

        y = (sin(x)/x-cos(x))/x**2

c-----------------------
      elseif(menu.eq.70) then
c-----------------------

        y = exp(-0.5*sin(x))

c-----------------------
      elseif(menu.eq.71) then
c-----------------------

        y = exp(-2.0*sin(x))

c-----------------------
      elseif(menu.eq.72) then
c-----------------------

        y = sin(x/2)**alpha72 * sin(p72*x)*sin(q72*x);

c-----------------------
      elseif(menu.eq.73) then
c-----------------------

        y = sin(x/2)**alpha73 * cos(pmq73*x);

c-----------------------
      elseif(menu.eq.80) then
c-----------------------

       srt = sqrt(3.0D0)
       sn = sin(x)
       cs = cos(x)
       f = srt/(sn+srt*cs)
       f4 = f**4
       f5 = f*f4
       y = (f5/5.0*sn-srt/4.0*f4+srt/5.0*f5*cs)*(sn-srt*cs)*sn

c------------
      end if
c-----------

c---
c done
c---

      return
      end
