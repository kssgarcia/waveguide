      subroutine integral_rec_fun (x,y,f,menu)

c========================================
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c-------------------------------
c compute (define) the integrand
c-------------------------------

      Implicit Double Precision (a-h,o-z)

      common/lgf/m,m1,m2

c-----------------------
      if(menu.eq.1) then
c-----------------------

        if((abs(x)+abs(y)).le.0.0000001) then
          x = 0.0001;
          y = 0.0001;
        end if
c       tmp = Dsin(0.5D0*x)**2 + Dsin(0.5D0*y)**2   ! square
        tmp = Dsin(0.5D0*x)**2   ! a38
     +      + Dsin(0.5D0*y)**2   ! a38
     +      + Dsin(0.5D0*(x+y))**2   ! a38
     +      + Dsin(0.5D0*(x-y))**2  ! a38 
        f = (1.0D0-Dcos(m1*x+m2*y))/tmp

c-----------------------
      elseif(menu.eq.2) then
c-----------------------

        if((abs(x)+abs(y)).le.0.0000001) then
          x = 0.001;
          y = 0.001;
        end if
        tmp = Dsin(0.5D0*x)**2 + Dsin(0.5D0*y)**2
        f = (1.0D0-Dcos(m1*x)*Dcos(m2*y))/tmp

c-----------------------
      elseif(menu.eq.3) then
c-----------------------

        if((abs(x)+abs(y)).le.0.0000001) then
          x = 0.001;
          y = 0.001;
        end if
        tmp = 1.0D0 - cos(x)*cos(y)
        if(abs(tmp).lt.0.0001) then
         tmp = 1.0;
        end if
        f = (1.0D0-cos(m1*x)*cos(m2*y))/tmp

c-----------------------
      elseif(menu.eq.4) then
c-----------------------

        tmp = 1.0D0-Dcos(1.0*x)*Dcos(1.0*y)
        if(abs(tmp).le.0.000001) then
         x = x+0.0001;
         y = y+0.0001;
         tmp = 1.0D0-Dcos(1.0D0*x)*Dcos(1.0D0*y)
        end if
        f = (1.0D0-Dcos(2.0D0*m*x))/tmp

c-----------------------
      elseif(menu.eq.5) then
c-----------------------

        if((abs(x)+abs(y)).le.0.0000001) then
          x = 0.001
          y = 0.001
        end if
        tmp = 1.0D0 - cos(x)*cos(y)
        if(abs(tmp).lt.0.0001) then
         tmp = 1.0
        end if
        arg = (m1+m2)*x+(m1-m2)*y
        f = (1.0D0-cos(arg))/tmp

c-----------------------
      elseif(menu.eq.10) then
c-----------------------

        if((abs(x)+abs(y)).le.0.0000001) then
          x = 0.0001D0
          y = 0.0001D0
        end if
        ari = 3-cos(m1*x+m2*y)-cos((m1-1)*x+m2*y)-cos(m1*x+(m2+1)*y)
        par = Dsin(0.5D0*x)**2   ! honey
     +      + Dsin(0.5D0*y)**2   !  honey
     +      + Dsin(0.5D0*(x+y))**2   !  honey
        f = ari/par

c-----------------------
      elseif(menu.eq.20) then
c-----------------------

        arg = 4.0D0-2.0D0*cos(x)-2.0D0*cos(y)

        if(abs(arg).le.0.0000001) then
          f = 0.0001D0
        else
          f = log(arg)
        end if

c-----------------------
      elseif(menu.eq.30) then
c-----------------------

        arg = 6.0D0-2.0D0*cos(x)-2.0D0*cos(y)-2.0D0*cos(x-y)

        if(abs(arg).le.0.0000001) then
          f = 0.0001D0
        else
          f = log(arg)
        end if

c-----------------------
      elseif(menu.eq.40) then
c-----------------------

        arg  = sin(x)**2+sin(y)**2
        arg1 = arg+1.0D0

        f = log(sqrt(arg1)+sqrt(arg))

c-----------------------
      elseif(menu.eq.41) then
c-----------------------

       f = x*12/32 + y*9/28 + x*y*(16-12-9)/(32*28)

c------------
      end if
c-----------

c---
c done
c---

      return
      end
