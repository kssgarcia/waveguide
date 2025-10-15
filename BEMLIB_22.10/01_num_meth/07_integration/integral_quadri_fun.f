      subroutine integral_quadri_fun (x,y,z,f,menu)

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c----------------------------------------

c-------------------------------
c compute (define) the integrand
c-------------------------------

      Implicit Double Precision (a-h,o-z)

      common/lgf/m,m1,m2,m3

c-----------------------
      if(menu.eq.1) then
c-----------------------

        if((abs(x)+abs(y)+abs(z)).le.0.0000001) then
          x = 0.001;
          y = 0.001;
          z = 0.001;
        end if
        tmp = Dsin(0.5*x)**2+sin(0.5*y)**2+sin(0.5*z)**2 ! bcc
     +       +Dsin(0.5*(x+y+z))**2
        tmp = Dsin(0.5*x)**2+sin(0.5*y)**2+sin(0.5*z)**2  ! fcc
     +       +Dsin(0.5*(x-y))**2
     +       +Dsin(0.5*(y-z))**2
     +       +Dsin(0.5*(z-x))**2
        f = (1.0D0-Dcos(m1*x+m2*y+m3*z))/tmp

c-----------------------
      elseif(menu.eq.2) then
c-----------------------

        if((abs(x)+abs(y)+abs(z)).le.0.0000001) then
          x = 0.001;
          y = 0.001;
          z = 0.001;
        end if
        tmp = Dsin(0.5*x)**2+sin(0.5*y)**2+sin(0.5*z)**2
     +       +Dsin(0.5*(x+y+z))**2
        f = (1.0D0-Dcos(m1*x)*Dcos(m2*y)*Dcos(m3*z))/tmp

c-----------------------
      elseif(menu.eq.3) then
c-----------------------

        if((abs(x)+abs(y)+abs(z)).le.0.0000001) then
          x = 0.001;
          y = 0.001;
          z = 0.001;
        end if
        tmp = 1.001D0-cos(x)*cos(y)*cos(z)   ! bcc
        tmp = 3.001D0-cos(x)*cos(y)-cos(y)*cos(z)-cos(z)*cos(x)  ! fcc
        f = (1.0D0-Dcos(m1*x+m2*y+m3*z))/tmp

c-----------------------
      elseif(menu.eq.4) then
c-----------------------

        if((abs(x)+abs(y)+abs(z)).le.0.0000001) then
          x = 0.001;
          y = 0.001;
          z = 0.001;
        end if
        tmp = 1.01D0-cos(x)*cos(y)*cos(z)
        tmp = 3.001D0-cos(x)*cos(y)-cos(y)*cos(z)-cos(z)*cos(x)  ! fcc
        f = (1.0D0-Dcos(m1*x)*Dcos(m2*y)*Dcos(m3*z))/tmp

c------------
      end if
c-----------

c---
c done
c---

      return
      end
