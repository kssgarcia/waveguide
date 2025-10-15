      subroutine davis 
     +
     +  (a,shrt
     +  ,x,sig,phi
     +  ,Ux,Uy,Uz
     +  )

c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licencing agreement.
c----------------------------------------

c-------------------------------
c Velocity for simple shear flow 
c over an orifice of zero thickness,
c using analytical expressions
c derived by AMJ Davis
c
c SYMBOLS:
c --------
c
c a:    orifice radius
c sig:  radial distance sigma in cylindrical coordinates
c phi:  azimuthal angle
c
c-------------------------------

      Implicit Double Precision (a-h,o-z)

c----------
c constants
c----------

      pi  = 3.14159265358
      pi3 = 3.0*pi

c---
c launching
c---

      sigs = sig**2
      phi2 = 2.0*phi

      as = a**2
      xs = x**2
      rs = sigs+xs

      rlms = 0.5*(rs-as+sqrt((rs-as)**2+4.0*xs*as))
      arlm = sqrt(rlms) 
      rlm  = arlm * x/abs(x)          ! lamda

      zet   = x/rlm       ! dimensionless
      ziets = zet**2      ! dimensionless
      xis   = rlms/as     ! dimensionless

      cs  = cos(phi)
      cs2 = cos(phi2)
      sn2 = sin(phi2)

      fc1 = a/pi3 * xis *zet * (1.0-zets)
     +                /((xis+1.0)*(xis+zets))

      fc2 = a/pi * zet * (1.0 - arlm/a * atan(a/arlm)
     +                        - zets/(3.0*(xis+zets)) 
     +                   )

      fc3 = 2.0*rlm/pi3 * zets/(xis+zets)
     +                  * sqrt((1.0-zets)/(1.0+xis))

      Uy = cs2*fc1 + fc2
      Uz = sn2*fc1
      Ux = cs *fc3

c------------------------------------------
c Add the simple shear flow above the plate
c------------------------------------------

      If(x.gt.0) Uy = Uy + shrt*x
     
c-----
c Done
c-----

  100 Format (1x,i3,3(1x,f15.10))

      Return
      End
