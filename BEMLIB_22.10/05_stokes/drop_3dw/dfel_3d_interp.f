      subroutine dfel_3d_interp
     +
     +    (vnx1,vny1,vnz1
     +    ,vnx2,vny2,vnz2
     +    ,vnx3,vny3,vnz3
     +    ,vnx4,vny4,vnz4
     +    ,vnx5,vny5,vnz5
     +    ,vnx6,vny6,vnz6
     +    ,al,be,ga
     +    ,xi,eta
     +
     +    ,vnx,vny,vnz
     +    )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c-----------------------------------
c Interpolate the unit normal vector
c-----------------------------------

      Implicit Double Precision (a-h,o-z)

c--------
c prepare
c--------

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

c----------------
c basis functions
c----------------

      ph2 = xi *(xi -al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi*(be+ga-1.0D0)/ga)/bec
      ph4 = xi *(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

c--------------------------
c compute the normal vector
c--------------------------

      vnx = vnx1*ph1 + vnx2*ph2 + vnx3*ph3 + vnx4*ph4
     +    + vnx5*ph5 + vnx6*ph6

      vny = vny1*ph1 + vny2*ph2 + vny3*ph3 + vny4*ph4
     +    + vny5*ph5 + vny6*ph6

      vnz = vnz1*ph1 + vnz2*ph2 + vnz3*ph3 + vnz4*ph4
     +    + vnz5*ph5 + vnz6*ph6

c-----
c Done
c-----

  100 Format (20(1x,f10.5))

      return
      end
