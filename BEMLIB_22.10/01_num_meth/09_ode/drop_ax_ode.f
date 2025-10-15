      subroutine drop_ax_ode
     +
     +  (npts
     +  ,capls
     +  ,Isp
     +  ,dpsi
     +  ,shp
     +  ,x,s
     +  ,volume
     +  )

c-----------------------------------------
c FDLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------

c--------------------------------------------
c Integrate ODEs by the modified Euler method
c with a uniform step size for the angle psi
c
c SYMBOLS:
c -------
c
c dpsi: increments in psi
c shp:  shooting parameter
c--------------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension x(256),s(256)

c--------
c prepare
c--------

      pi  = 3.14159 265358D0
      dpsih = 0.5D0*dpsi

c----------------
c top of the drop
c----------------

      psi  = 0.0D0
      x(1) = 0.0D0
      s(1) = 0.0D0

      Do i=1,npts

       If(i.eq.1) then
        xp = 0.0D0              ! from eqs. (4.4.17)
        sp = 2.0D0/shp          ! from eqs. (4.4.17)
       Else
        Q  = Dsin(psi)/s(i)+Isp*x(i)/capls-shp
        xp = Dsin(psi)/Q
        sp =-Dcos(psi)/Q
       End If

       xsv  = x(i)        ! save
       ssv  = s(i)
       xpsv = xp
       spsv = sp

       i1 = i+1

       x(i1) = x(i)+xp*dpsi
       s(i1) = s(i)+sp*dpsi

       psi = psi+dpsi

       Q  = Dsin(psi)/s(i1)+Isp*x(i1)/capls-shp
       xp = Dsin(psi)/Q
       sp =-Dcos(psi)/Q

       x(i1) = xsv + (xpsv+xp)*dpsih
       s(i1) = ssv + (spsv+sp)*dpsih

      End Do

c-------------------------------------------
c compute the volume of the integrated shape
c by the trapezoidal rule
c-------------------------------------------

      volume = 0.0D0

      Do i=1,npts
       volume = volume+(s(i+1)**2+s(i)**2)*Dabs(x(i+1)-x(i))
      End Do

      volume = 0.5D0*volume     ! to account for trapezoidal weights
      volume = pi*volume

c-----
c done
c-----

      Return
      End
